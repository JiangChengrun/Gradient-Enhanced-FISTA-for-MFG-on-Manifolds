% This function tests the MF-Fista algorithm for the 2D membrane free-energy functional






close all;
clear;

% sphere in 3d
[surf.pt,surf.trg] = readOFF(['../Data/sphere.off']);
surf = surfOperators(surf);

% sphere in 5d
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
opts.funcF = @(rho) zeros(size(rho));
opts.gradF = @(rho) zeros(size(rho));
lambdaG = 8e-1;
opts.funcG = @(rhoend) lambdaG*rhoend.*(log(rhoend)-logrho1);
opts.gradG = @(rhoend) lambdaG*(1+log(rhoend)-logrho1);

opts.plot = 1;
opts.savegif = 0;
opts.saveshot = 0;
opts.nt = 32;

opts.maxit = 3000;
opts.tol = 1e-5;

opts.stepsize0 = 8e0;
opts.stepmodif = 0.5;
opts.submaxit = 5;
opts.acc = 0;

%% FISTA

tic;
[rho, flux, output] = mfgMfFista(surf,rho0,opts);
toc;

tic;
[rho5d, flux5d, output5d] = mfgMfFista(surf5d,rho0,opts);
toc;

%% compare
filenameSave = 'sphere_compare';
save(['results/',filenameSave]);

% cost
obj = output.objArray(end);
obj5d = output5d.objArray(end);
fprintf('difference between obj: %e \n',abs(obj(end)-obj5d));

figure;hold on
plot(output.objArray,'r--','LineWidth',2,'DisplayName','3D');
plot(output5d.objArray,'b:','LineWidth',2,'DisplayName','5D');
legend()
set(gcf,'unit','centimeters','position',[10 5 7 7])
set(gca,'Position',[0.2,0.2,0.7,0.7]);
xlabel('number of iterations');ylabel('total cost');
% title('objective function value');
fig = gcf;
exportgraphics(fig,['results\',filenameSave,'_obj.eps']);

% optimizer
flux_3d25d = cat(3,flux,zeros(nTrg,opts.nt,2));
for t = 1:opts.nt
    flux_t = squeeze(flux_3d25d(:,t,:)); 
    flux_t5d = flux_t*U;
    flux_3d25d(:,t,:) = reshape(flux_t5d,nTrg,1,5);
end
fprintf('difference between density: %e \n',max(abs(rho(:)-rho5d(:))));
fprintf('difference between flux: %e \n',max(abs(flux_3d25d(:)-flux5d(:))));


