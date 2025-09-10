% test_eight
% This function tests the MF-Fista algorithm for the 2D membrane free-energy functional
% Input: None
% Output: None







close all;
clear;

meshName = 'eight';
load(['../Data/',meshName,'.mat',])
surf.pt = pt(:,[3,1,2]);
surf.trg =  trg;
clear pt trg

nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
surf = surfOperators(surf);

%% example
% distance = vecnorm(surf.pt,2,2);
rho0 = zeros(nPt,1);
rho1 = zeros(nPt,1);

distance = vecnorm(surf.pt-surf.pt(522,:),2,2);
rho0 = 25*exp(-distance.^2/0.01);
distance = vecnorm(surf.pt-surf.pt(330,:),2,2);
rho1 = 25*exp(-distance.^2/0.01);

% aroundpt = getaroundpt(surf.pt,surf.trg);
% rho0(aroundpt{522}) = 1e2;
% rho1(aroundpt{330}) = 1e2; 

% rho0 = rho0 + 0.1;
% rho1 = rho1 + 0.1;
% logrho1 = log(rho1);

figure(1);clf;viewMesh(surf,rho0);axis on; grid on; view(3);colorbar
figure(2);clf;viewMesh(surf,rho1);axis on; grid on; view(3);colorbar

% x = surf.trgCenter(:,1); y = surf.trgCenter(:,2); z = surf.trgCenter(:,3);
% obstacle = zeros(nTrg,1);
x = surf.pt(:,1); y = surf.pt(:,2); z = surf.pt(:,3);
obstacle = zeros(nPt,1);
obstacle(x> 0.19 & x< 0.31 & y>0) = 1;
obstacle(x<-0.19 & x>-0.31 & y<0) = 1;
% rho0(obstacle>0) = 0;
% rho1(obstacle>0) = 0;
% figure(3);clf;viewMesh(surf,obstacle); axis on; grid on; view(3)

egName = 'obs';
obsPenalty = 50;
opts.funcL = @(rho,m) sum(m.^2,3)./(2*rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
opts.gradLm = @(rho,m) m./rho; 
opts.funcF = @(rho) obsPenalty.*obstacle.*rho;
opts.gradF = @(rho) obsPenalty.*obstacle;
lambdaG = 10;
opts.funcG = @(rhoend) lambdaG/2*(rhoend-rho1).^2;
opts.gradG = @(rhoend) lambdaG*(rhoend-rho1);

% parameters
opts.plot = 1;
opts.savegif = 1;
opts.saveshot = 1;

opts.nt = 32;
opts.maxit = 5e3;
opts.tol = 1e-7;

opts.stepsize0 = 1e1;
opts.stepmodif = 0.5;
opts.submaxit = 5;
opts.acc = 1;

%% FISTA

tic;
% [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
[rho, flux, output] = mfgMfFista(surf,rho0,opts);
totaltime = toc
totalnit = length(output.objArray)
ittime = totaltime/totalnit

rhomin = min(rho,[],'all');
rhomax = max(rho,[],'all');

rhos = surf.pt2trg*(rho(:,1:end-1)+rho(:,2:end))/2;
cost = output.costArray{end};
fprintf('dynamic cost: %f \n',...
        sum(surf.trgArea.*opts.funcL(rhos,flux),'all')/opts.nt);
fprintf('interaction cost: %f \n',...
        sum(surf.ptArea.*opts.funcF(rho(:,2:end-1)),'all')/opts.nt);
fprintf('terminal cost: %f \n',...
        sum(surf.ptArea.*opts.funcG(rho(:,end)),'all'));
fprintf('Total cost: %f \n',output.objArray(end));
figure;plot(sum(rho.*surf.ptArea));


% visualization
filenameSave = [meshName,'_mfg_',egName];
save(['results/',filenameSave]);

%%
isObs = min(ismember(surf.trg,find(obstacle>0)),[],2);

close all
if opts.plot
    %-------show and save objective history
    %h = figure(4);
    figure(1);
    plot(output.objArray,'LineWidth',2);
    print('-dpng',['results/' filenameSave,'_obj.png']);
    
    %-------show evolution
    figure(2);
    set(gcf,'color','w');
    for i = 1:size(rho,2)-1
        clf
        
        patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
        set(patchNonobs,'FaceVertexCData',rho(:,i));
        set(patchNonobs,'FaceColor','interp');
        set(patchNonobs,'edgecolor','none');
        
        patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
        set(patchObs,'FaceVertexCData',rho(:,i));
        set(patchObs,'FaceColor','interp');
        set(patchObs,'FaceAlpha',0.3)
        set(patchObs,'edgecolor','none');

        caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
%         viewVectF(surf.trgCenter,squeeze(flux(:,i,:)));

        daspect([1 1 1])
        view(3)
        axis tight
        axis off
        axis equal
        material dull
        camlight
        lighting phong
        
        frame = getframe(gcf);
        im{i} = frame2im(frame);
        pause(0.2);        
    end
    i = size(rho,2);
    clf

    patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
    set(patchNonobs,'FaceVertexCData',rho(:,i));
    set(patchNonobs,'FaceColor','interp');
    set(patchNonobs,'edgecolor','none');
    
    patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
    set(patchObs,'FaceVertexCData',rho(:,i));
    set(patchObs,'FaceColor','interp');
    set(patchObs,'FaceAlpha',0.3)
    set(patchObs,'edgecolor','none');

    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on

    daspect([1 1 1])
    view(3)
    axis tight
    axis off
    axis equal
    material dull
    camlight
    lighting phong
    
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
        if ~isfield(opts,'num_frame') num_frame = 10; end
        idx_frame = round(linspace(1,opts.nt+1,num_frame));
        t_frame = (idx_frame-1)./(opts.nt);

        close all
        for k = 1:num_frame-1
            clf
            idx = idx_frame(k);

            patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
            set(patchNonobs,'FaceVertexCData',rho(:,idx));
            set(patchNonobs,'FaceColor','interp');
            set(patchNonobs,'edgecolor','none');

            patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
            set(patchObs,'FaceVertexCData',rho(:,idx));
            set(patchObs,'FaceColor','interp');
            set(patchObs,'FaceAlpha',0.3)
            set(patchObs,'edgecolor','none');
            
            caxis([min(rho(:,idx)),max(rho(:,idx))]);colorbar;hold on
%             viewVectF(surf.trgCenter,squeeze(flux(:,idx,:)));

            daspect([1 1 1])
            view(3)
            axis tight
            axis off
            axis equal
            material dull
            camlight
            lighting phong

            set(gcf,'unit','centimeters','position',[10 5 2 3])
            set(gca,'Position',[0.1,0.05,0.65,0.8]);
            caxis([min(rho(:,idx)),max(rho(:,idx))]);
            colorbar('Position',[0.76,0.1,0.05,0.7]);
            title(['t=',num2str(t_frame(k))]);
            fig = gcf;
            exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        end
        clf
        k = num_frame;
        idx = idx_frame(k);

        patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
        set(patchNonobs,'FaceVertexCData',rho(:,idx));
        set(patchNonobs,'FaceColor','interp');
        set(patchNonobs,'edgecolor','none');

        patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
        set(patchObs,'FaceVertexCData',rho(:,idx));
        set(patchObs,'FaceColor','interp');
        set(patchObs,'FaceAlpha',0.3)
        set(patchObs,'edgecolor','none');        
        
        caxis([min(rho(:,idx)),max(rho(:,idx))]);colorbar;hold on

        daspect([1 1 1])
        view(3)
        axis tight
        axis off
        axis equal
        material dull
        camlight
        lighting phong
        
        set(gcf,'unit','centimeters','position',[10 5 2 3])
        set(gca,'Position',[0.1,0.05,0.65,0.8]);
        caxis([min(rho(:,idx)),max(rho(:,idx))]);
        colorbar('Position',[0.76,0.1,0.05,0.7]);
        title(['t=',num2str(t_frame(k))]);
        fig = gcf;
        exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        
        % different view angle
        for k = 1:num_frame-1
            clf
            idx = idx_frame(k);
            
            patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
            set(patchNonobs,'FaceVertexCData',rho(:,idx));
            set(patchNonobs,'FaceColor','interp');
            set(patchNonobs,'edgecolor','none');

            patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
            set(patchObs,'FaceVertexCData',rho(:,idx));
            set(patchObs,'FaceColor','interp');
            set(patchObs,'FaceAlpha',0.3)
            set(patchObs,'edgecolor','none');        
            
            caxis([min(rho(:,idx)),max(rho(:,idx))]);colorbar;hold on
%             viewVectF(surf.trgCenter,squeeze(flux(:,idx,:)));

            daspect([1 1 1])
            view(-90,90)
            axis tight
            axis off
            axis equal
            material dull
            camlight
            lighting phong

            set(gcf,'unit','centimeters','position',[10 5 2 3])
            set(gca,'Position',[0.1,0.05,0.65,0.8]);
            caxis([min(rho(:,idx)),max(rho(:,idx))]);
            colorbar('Position',[0.76,0.1,0.05,0.7]);
            title(['t=',num2str(t_frame(k))]);
            fig = gcf;
            exportgraphics(fig,['results/',filenameSave,'_a2_shot',num2str(k),'.eps']);
        end
        clf
        k = num_frame;
        idx = idx_frame(k);
        
        patchNonobs = patch('Vertices',surf.pt, 'Faces',surf.trg(~isObs,:));
        set(patchNonobs,'FaceVertexCData',rho(:,idx));
        set(patchNonobs,'FaceColor','interp');
        set(patchNonobs,'edgecolor','none');

        patchObs = patch('Vertices',surf.pt, 'Faces',surf.trg(isObs,:));
        set(patchObs,'FaceVertexCData',rho(:,idx));
        set(patchObs,'FaceColor','interp');
        set(patchObs,'FaceAlpha',0.3)
        set(patchObs,'edgecolor','none');        

        caxis([min(rho(:,idx)),max(rho(:,idx))]);colorbar;hold on

        daspect([1 1 1])
        view(-90,90)
        axis tight
        axis off
        axis equal
        material dull
        camlight
        lighting phong
        
        set(gcf,'unit','centimeters','position',[10 5 2 3])
        set(gca,'Position',[0.1,0.05,0.65,0.8]);
        caxis([min(rho(:,idx)),max(rho(:,idx))]);
        colorbar('Position',[0.76,0.1,0.05,0.7]);
        title(['t=',num2str(t_frame(k))]);
        fig = gcf;
        exportgraphics(fig,['results/',filenameSave,'_a2_shot',num2str(k),'.eps']);
    
    end
    
end
