% This function tests the MF-Fista algorithm for the 2D membrane free-energy functional

close all;
clear;

meshName = 'mobius_reg64';
load(['../Data/',meshName,'.mat'],'surf');

nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
surf = surfOperators(surf);

%% example
% aroundpt = getaroundpt(surf.pt,surf.trg);
% rho0 = exp(-(vecnorm(surf.pt-[0.8,0.8,0],2,2)).^2/0.01);
% rho1 = exp(-(vecnorm(surf.pt-[0.2,0.2,0],2,2)).^2/0.01);
% 
% rho0 = rho0./(sum(rho0.*surf.ptArea)) + 0.1;
% 
% 
% figure;viewMesh(surf,rho0);colorbar
% figure;viewMesh(surf,rho1);colorbar

sigma = 0.3;
rho0 = 10*exp(-(vecnorm(surf.flatpt-[2*pi/4,0],2,2)).^2/0.04);
rho1 = 10*exp(-(vecnorm(surf.flatpt-[6*pi/4,-0],2,2)).^2/0.04);
% rho0 = zeros(nPt,1);
% rho1 = zeros(nPt,1);
% rho0(surf.flatpt(:,1)>1*pi/4-0.2 & surf.flatpt(:,1)<1*pi/4+0.2) = 10;
% rho1(surf.flatpt(:,1)>7*pi/4-0.2 & surf.flatpt(:,1)<7*pi/4+0.2) = 10;
rho1 = rho1*(sum(rho0.*surf.ptArea)/sum(rho1.*surf.ptArea));

rho0 = rho0 + 0.1;
rho1 = rho1 + 0.1;
logrho1 = log(rho1);

figure(1);viewMesh(surf,rho0);colorbar
figure(2);viewMesh(surf,rho1);colorbar


%% parameters
opts.plot = 1;
opts.savegif = 1;
opts.saveshot = 0;
opts.funcL = @(rho,m) sum(m.^2,3)./(2*rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
opts.gradLm = @(rho,m) m./rho; 
lambdaF = 2;
% egName = 'vanilla2';
% opts.funcF = @(rho) zeros(size(rho));
% opts.gradF = @(rho) zeros(size(rho));
egName = 'rhologrho';
opts.funcF = @(rho) lambdaF*rho.*log(rho);
opts.gradF = @(rho) lambdaF*(log(rho)+1);
lambdaG = 1e2;
opts.funcG = @(rhoend) lambdaG*rhoend.*(log(rhoend)-logrho1);
opts.gradG = @(rhoend) lambdaG*(log(rhoend)-logrho1+1);

opts.nt = 32;

opts.maxit = 40000;
opts.tol = 1e-6;

opts.stepsize0 = 1e-1;
opts.stepmodif = 0.5;
opts.submaxit = 5;

%% FISTA

tic;
[rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
% [rho, flux, output] = mfgMfFista(surf,rho0,opts);
totaltime = toc
totalnit = length(output.objArray)
ittime = totaltime/totalnit

%
rhomin = min(rho,[],'all');
rhomax = max(rho,[],'all');

rhos = surf.pt2trg*(rho(:,1:end-1)+rho(:,2:end))/2;
% cost = output.costArray{end};

%%
rhomin = min(rho,[],'all')-0.1;
rhomax = max(rho,[],'all')+0.1;
fprintf('Wasserstein 2 distance: %f \n',output.objArray(end));

fprintf('dynamic cost: %f \n',...
        sum(surf.trgArea.*opts.funcL(rhos,flux),'all')/opts.nt);
fprintf('interaction cost: %f \n',...
        sum(surf.ptArea.*opts.funcF(rho(:,2:end-1)),'all')/opts.nt);
fprintf('terminal cost: %f \n',...
        sum(surf.ptArea.*opts.funcG(rho(:,end)),'all'));
fprintf('Total cost: %f \n',output.objArray(end));
figure;plot(sum(rho.*surf.ptArea));

filenameSave = [meshName,'_mfg_',egName];
save(['results/',filenameSave]);

%% visualization
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
        viewMesh(surf,rho(:,i));view(50,50);
        caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
        frame = getframe(gcf);
        im{i} = frame2im(frame);
        pause(0.2);        
    end
    i = size(rho,2);
    clf
    viewMesh(surf,rho(:,i));view(50,50);
    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
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
%
    %---------- save eps
    if opts.saveshot
        if ~isfield(opts,'num_frame') num_frame = 8; end
        idx_frame = round(linspace(1,opts.nt+1,num_frame));
        t_frame = (idx_frame-1)./(opts.nt);

        close all
        for k = 1:num_frame-1
            clf
            idx = idx_frame(k);
%             fluxlength = vecnorm(flux(:,idx,:),2,3);
%             [~,showvecidx] = maxk(fluxlength,5);
%             [~,sortvecidx] = sort(fluxlength);
            viewMesh(surf,rho(:,idx));view(50,50);hold on;%view(10,40);
%             viewVectF(surf.trgCenter(showvecidx,:),squeeze(flux(showvecidx,idx,:)));
%             viewVectF(surf.trgCenter(sortvecidx(1:20:end),:),squeeze(flux(sortvecidx(1:20:end),idx,:)));
            set(gcf,'unit','centimeters','position',[10 5 2 3])
            set(gca,'Position',[0.1,0.05,0.65,0.8]);% left margin, lower margin, width, height
            caxis([min(rho(:,idx)),max(rho(:,idx))]);
            colorbar('Position',[0.76,0.1,0.05,0.7]);% left margin, lower margin, width,
            title(['t=',num2str(t_frame(k))]);
            fig = gcf;
            exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        end
        clf
        k = num_frame;
        idx = idx_frame(k);
        viewMesh(surf,rho(:,idx));view(50,50);hold on;%view(10,40);
        set(gcf,'unit','centimeters','position',[10 5 2 3])
        set(gca,'Position',[0.1,0.05,0.65,0.8]);
        caxis([min(rho(:,idx)),max(rho(:,idx))]);
        colorbar('Position',[0.76,0.1,0.05,0.7]);
        title(['t=',num2str(t_frame(k))]);
        fig = gcf;
        exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        
    end
    
end


