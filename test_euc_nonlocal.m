% This function tests the MF-Fista algorithm for the 2D membrane free-energy functional
% Input: None
% Output: None
close all;
clear;

meshName = 'euc_equ_reg24';
% meshName = 'euc_right_reg24';
load(['../Data/',meshName,'.mat'],'surf');

nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
surf = surfOperators(surf);

%% example
aroundpt = getaroundpt(surf.pt,surf.trg);
rho0 = exp(-(vecnorm(surf.pt-[0.9,0.5,0],2,2)).^2/0.01);
rho1 = 2*exp(-10*(surf.pt(:,2)-0.5).^2-(surf.pt(:,1)-0.1).^2).*((surf.pt(:,1)-0.1).^2-1);

rho0 = rho0./(sum(rho0.*surf.ptArea)) + 0.1;


figure;viewMesh(surf,rho0);colorbar
figure;viewMesh(surf,rho1);colorbar

%% kernel
mu = 1;
sigmasq = 0.01;
kernel = sparse(nPt,nPt);
for iPt = 1:nPt
    firstringPt = aroundpt{iPt};
    kernel(iPt,firstringPt) = mu*exp(-vecnorm( surf.pt(firstringPt,:)-surf.pt(iPt,:),2,2 ).^2/sigmasq);
end

%% parameters
klPenalty = 3;
opts.plot = 1;
opts.savegif = 1;
egName = 'nonlocal';
opts.funcL = @(rho,m) sum(m.^2,3)./(2*rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
opts.gradLm = @(rho,m) m./rho; 
opts.funcF = @(rho) 1/2*trace(rho'*kernel*rho);
opts.gradF = @(rho) 1*kernel*rho;
opts.funcG = @(rhoend) 5*rhoend.*rho1;
opts.gradG = @(rhoend) 5*rho1;

opts.nt = 32;

opts.maxit = 3000;
opts.tol = 1e-3;

opts.stepsize0 = 10;
opts.stepmodif = 0.5;
opts.submaxit = 5;



%% FISTA

tic;
% [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
[rho, flux, output] = mfgMfFista(surf,rho0,opts);
toc;

%%
rhomin = min(rho,[],'all')-0.1;
rhomax = max(rho,[],'all')+0.1;
fprintf('Wasserstein 2 distance: %f \n',output.objArray(end));
figure;plot(sum(rho.*surf.ptArea));


% %% visualization
% filenameSave = ['mfg_',meshName,'_',num2str(opts.nt),egName,'_',num2str(opts.maxit)];
% 
% if opts.plot
%     %h = figure(4);
%     figure(1);
%     plot(output.objArray,'LineWidth',2);
%     print('-dpng',['results/' filenameSave,'_obj.png']);
%     figure(2);
% %     ha = tight_subplot(1,2,[.01 .03],[.1 .01],[.01 .01]);
%     set(gcf,'color','w');
% %     set(gca,'position',[0 0 1 1],'units','normalized')
%     for i = 1:size(rho,2)
%         trimesh(surf.trg,surf.pt(:,1),surf.pt(:,2),rho(:,i));
%         axis([0 1 0 1 rhomin rhomax])
%         frame = getframe(gcf);
%         im{i} = frame2im(frame);
%         pause(0.2);        
%     end
%     
%    if opts.savegif
%         
%         for idx = 1:size(rho,2)
%             [A,map] = rgb2ind(im{idx},256);
%             if idx == 1
%                 imwrite(A,map,['results/',filenameSave,'.gif'],...
%                     'gif','LoopCount',Inf,'DelayTime',0.2);
%             else
%                 imwrite(A,map,['results/',filenameSave,'.gif'],...
%                     'gif','WriteMode','append','DelayTime',0.2);
%             end
%         end
%     end
% 
% end
% save(['results/',filenameSave,]);

%% visualization
filenameSave = [meshName,'_',num2str(opts.nt),egName,'_',num2str(opts.maxit)];

im = cell(size(rho,2),1);
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
        pause(0.2);        
    end
    i = size(rho,2);
    clf;
    viewMesh(surf,rho(:,i));
    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
    frame = getframe(gcf);
    im{i} = frame2im(frame);
    
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

