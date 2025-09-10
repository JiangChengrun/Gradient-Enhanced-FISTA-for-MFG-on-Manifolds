% test_kitten
% This function tests the MF-Fista algorithm for the 2D membrane free-energy functional






close all;
clear;

meshName = 'genus3';
[surf.pt,~,surf.trg] = readObjShape(['../../Data/',meshName,'.obj']);

nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
surf = surfOperators(surf);

% load('results/genus3_32_congest_6000','-mat','rho','flux');
% opts.rho = rho(:,2:end);
% opts.m = flux;

%% example
% distance = vecnorm(surf.pt,2,2);
x = surf.pt(:,1); y = surf.pt(:,2); z = surf.pt(:,3);
rho0 = zeros(nPt,1);
rho1 = zeros(nPt,1);

distance = vecnorm(surf.pt-surf.pt(2057,:),2,2);
rho0 = 50*exp(-distance.^2).*(distance<0.07);
distance = vecnorm(surf.pt-surf.pt(5753,:),2,2);
rho1 = 50*exp(-distance.^2).*(distance<0.07);

% aroundpt = getaroundpt(surf.pt,surf.trg);
% rho0(aroundpt{1758}) = 1e3;
% rho1(aroundpt{2057}) = 1e3; %6061
% rho1(aroundpt{1118}) = 1e3;
% rho1(aroundpt{3169}) = 1e3;

rho0 = rho0 + 0.1;
rho1 = rho1 + 0.1;
logrho1 = log(rho1);

figure(1);clf;viewMesh(surf,rho0);axis on; grid on; view(2)
figure(2);clf;viewMesh(surf,rho1);axis on; grid on; view(2)

x = surf.trgCenter(:,1); y = surf.trgCenter(:,2); z = surf.trgCenter(:,3);
obstacle = zeros(nTrg,1);
obstacle(y>0&y<0.2&z>0) = 1;
figure(3);clf;viewMesh(surf,obstacle); axis on; grid on; view(2)

egName = '_obs';
obsPenalty = 5e4;
aggrePenalty = 2;
opts.funcF = @(rho,m) sum(m.^2,3)./(2*rho) + obsPenalty.*obstacle.*rho;
opts.funcGrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2) + obsPenalty.*obstacle;
opts.funcGm   = @(rho,m) m./rho;
% egName = '';
% opts.funcF = @(rho,m) sum(m.^2,3)./(2*rho);
% opts.funcGrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
% opts.funcGm   = @(rho,m) m./rho;
% opts.funcFend = @(rhoend) klPenalty* rhoend.*(log(rhoend)-logrho1).*(rhoend > 0);
% opts.funcGrhoend = @(rhoend) klPenalty* (log(rhoend)-logrho1 + 1).*(rhoend > 0);
klPenalty = 1;
opts.funcFend = @(rhoend) klPenalty/2* (rhoend-rho1).^2;
opts.funcGrhoend = @(rhoend) klPenalty* (rhoend-rho1);

%% parameters
opts.plot = 1;
opts.savegif = 1;
opts.nt = 32;

opts.maxit = 2000;
opts.tol = 1e-3;

opts.stepsize0 = 0.01;
opts.stepmodif = 0.5;
opts.submaxit = 5;



%% FISTA

tic;
% [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
[rho, flux, output] = mfgMfFista(surf,rho0,opts);
toc;

%%
zmin = min(rho,[],'all')-0.1;
zmax = max(rho,[],'all')+0.1;
fprintf('Wasserstein 2 distance: %f \n',output.objArray(end));


%% visualization
filenameSave = [meshName,'_',num2str(opts.nt),egName,'_',num2str(opts.maxit)];

if opts.plot
    %h = figure(4);
    figure(1);
    plot(output.objArray,'LineWidth',2);
    print('-dpng',['results/' filenameSave,'_obj.png']);

    figure(2);
%     ha = tight_subplot(1,2,[.01 .03],[.1 .01],[.01 .01]);
    set(gcf,'color','w');
%     set(gca,'position',[0 0 1 1],'units','normalized')
    for i = 1:size(rho,2)-1
        clf
        viewMesh(surf,rho(:,i));
        caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
        viewVectF(surf.trgCenter,squeeze(flux(:,i,:)));
        frame = getframe(gcf);
        im{i} = frame2im(frame);
        pause(0.5);        
    end
    i = size(rho,2);
    clf;
    viewMesh(surf,rho(:,i));
    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
    frame = getframe(gcf);
    im{i} = frame2im(frame);
    
    figure(3);plot(sum(rho.*surf.ptArea));
    
   if opts.savegif
        
        for idx = 1:size(rho,2)
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,['results/',filenameSave,'.gif'],...
                    'gif','LoopCount',Inf,'DelayTime',0.5);
            else
                imwrite(A,map,['results/',filenameSave,'.gif'],...
                    'gif','WriteMode','append','DelayTime',0.5);
            end
        end
    end

end
save(['results/',filenameSave,]);

