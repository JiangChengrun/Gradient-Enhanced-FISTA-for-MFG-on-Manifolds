% test_nonlocal

%function 
close all;
clear;

% sphere in 3d
[surf.pt,surf.trg] = readOFF(['../Data/sphere.off']);

% sphere in 5d
meshName = 'sphere5d';
rng('default')
A = rand(5,5);
[U,S,V] = svd(A);
nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
surf5d.pt = [surf.pt,zeros(nPt,2)]*U;
surf5d.trg = surf.trg;
surf5d = surfOperators5d(surf5d);

%% example
aroundpt = getaroundpt(surf.pt,surf.trg);
% rho0 = zeros(nPt,1);
% rho1 = zeros(nPt,1);
rho0 = exp(-(vecnorm(surf.pt-[cos(pi/3),0,sin(pi/3)],2,2)).^2/0.05);
rho1 = exp(-(vecnorm(surf.pt-[cos(2*pi/3),0,sin(2*pi/3)],2,2)).^2/0.05);

rho0 = rho0 + 0.1;
rho1 = rho1 + 0.1;
logrho1 = log(rho1);


figure(1);viewMesh(surf,rho0);colorbar
figure(2);viewMesh(surf,rho1);colorbar


%% parameters

opts.funcL = @(rho,m) sum(m.^2,3)./(2*rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
opts.gradLm = @(rho,m) m./rho; 
egName = 'vanilla';
opts.funcF = @(rho) zeros(size(rho));
opts.gradF = @(rho) zeros(size(rho));
% egName = 'sqrtrho';
% lambdaF = 2e-1;
% opts.funcF = @(rho) lambdaF*sqrt(rho+1e-4);
% opts.gradF = @(rho) lambdaF/2./sqrt(rho+1e-4);
lambdaG = 8e-1;
opts.funcG = @(rhoend) lambdaG*rhoend.*(log(rhoend)-logrho1);
opts.gradG = @(rhoend) lambdaG*(1+log(rhoend)-logrho1);


opts.plot = 1;
opts.savegif = 1;
opts.saveshot = 1;
opts.nt = 32;

opts.maxit = 3000;
opts.tol = 1e-5;

opts.stepsize0 = 8e0;
opts.stepmodif = 0.5;
opts.submaxit = 5;
opts.acc = 0;

%% FISTA

tic;
% [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
[rho, flux, output] = mfgMfFista(surf5d,rho0,opts);
toc;

%
rhomin = min(rho,[],'all')-0.1;
rhomax = max(rho,[],'all')+0.1;

rhos = surf5d.pt2trg*(rho(:,1:end-1)+rho(:,2:end))/2;
cost = output.costArray{end};

fprintf('dynamic cost: %f \n',...
        sum(surf5d.trgArea.*opts.funcL(rhos,flux),'all')/opts.nt);
fprintf('interaction cost: %f \n',...
        sum(surf5d.ptArea.*opts.funcF(rho(:,2:end-1)),'all')/opts.nt);
fprintf('terminal cost: %f \n',...
        sum(surf5d.ptArea.*opts.funcG(rho(:,end)),'all'));
fprintf('Total cost: %f \n',output.objArray(end));
figure;plot(sum(rho.*surf5d.ptArea));



%% visualization
filenameSave = [meshName,'_mfg_',egName];
save(['results/',filenameSave]);

close all
if opts.plot
    %-------show and save objective history
    %h = figure(4);
    figure(1);
    plot(output.objArray,'LineWidth',2);
    print('-dpng',['results/' filenameSave,'_obj.png']);
    
    %-------show evolution
    figure(2);
%     ha = tight_subplot(1,2,[.01 .03],[.1 .01],[.01 .01]);
    set(gcf,'color','w');
%     set(gca,'position',[0 0 1 1],'units','normalized')
    for i = 1:size(rho,2)-1
        clf
        viewMesh(surf,rho(:,i));
        caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
%         caxis([zmin,zmax]);colorbar;hold on
%         viewVectF(surf.trgCenter,squeeze(flux(:,i,:)));
        frame = getframe(gcf);
        im{i} = frame2im(frame);
        pause(0.2);        
    end
    i = size(rho,2);
    clf
    viewMesh(surf,rho(:,i));
    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
%     caxis([zmin,zmax]);colorbar;hold on
    frame = getframe(gcf);
    im{i} = frame2im(frame);
    pause(0.2);
    
    %--------- save gif
    if opts.savegif
        for idx = 1:size(rho,2)
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,['results/',filenameSave,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.3);
            else
                imwrite(A,map,['results/',filenameSave,'.gif'],'gif','WriteMode','append','DelayTime',0.3);
            end
        end
    end

    %---------- save eps
    if opts.saveshot
        if ~isfield(opts,'num_frame') num_frame = 5; end
        idx_frame = round(linspace(1,opts.nt+1,num_frame));
        t_frame = (idx_frame-1)./(opts.nt);

        close all
        for k = 1:num_frame-1
            clf
            idx = idx_frame(k);
            viewMesh(surf,rho(:,idx));hold on;
%             viewVectF(surf.trgCenter(1:5:end,:),squeeze(flux(1:5:end,idx,:)));
            set(gcf,'unit','centimeters','position',[10 5 2 3])
            set(gca,'Position',[0.1,0.05,0.65,0.8]);% left margin, lower margin, width, height
            caxis([min(rho(:,idx)),max(rho(:,idx))]);
%             caxis([zmin,zmax]);colorbar;hold on
            colorbar('Position',[0.76,0.1,0.05,0.7]);% left margin, lower margin, width,
            title(['t=',num2str(t_frame(k))]);
            fig = gcf;
            exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        end
        clf
        k = num_frame;
        idx = idx_frame(k);
        viewMesh(surf,rho(:,idx));hold on;
%         viewVectF(surf.trgCenter(1:5:end,:),squeeze(flux(1:5:end,idx,:)));
        set(gcf,'unit','centimeters','position',[10 5 2 3])
        set(gca,'Position',[0.1,0.05,0.65,0.8]);% left margin, lower margin, width, height
        caxis([min(rho(:,idx)),max(rho(:,idx))]);
%         caxis([zmin,zmax]);colorbar;hold on
        colorbar('Position',[0.76,0.1,0.05,0.7]);% left margin, lower margin, width,
        title(['t=',num2str(t_frame(k))]);
        fig = gcf;
        exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        
    end
    
end


