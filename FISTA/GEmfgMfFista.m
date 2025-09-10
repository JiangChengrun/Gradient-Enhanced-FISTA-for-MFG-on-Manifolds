function [rho,m,outs] = GEmfgMfFista(surf,rho0,opts)
    % This function computes the MF-Fista algorithm for minimizing the
    % 2D membrane free-energy functional
    % Input: surf = surface operator
    %        rho0 = initial value of rho
    %        opts = options
    % Output: rho = final value of rho
    %         m = final value of m
    %         outs = output structure   
%% parameters
nPt = size(surf.pt,1); %number of points
nTrg = size(surf.trg,1); %number of triangles
nDim = size(surf.pt,2); %number of dimension
% idIso = find(surf.ptArea==0);

if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
dt = 1/nt;


if isfield(opts,'acc') acc = logical(opts.acc); else acc = logical(0); end
% inital value of rho and m
if isfield(opts,'rho') xrhoNitm = opts.rho; else xrhoNitm = ones(nPt,nt); end
if isfield(opts,'m')   xmNitm  = opts.m;    else xmNitm  = zeros(nTrg,nt,nDim); end

% functions and gradients
if isfield(opts,'funcL') funcL = opts.funcL; else funcL = @(rho,m) sum(m.^2,3)./(2*rho); end
if isfield(opts,'gradLrho') gradLrho = opts.gradLrho; else gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2); end
if isfield(opts,'gradLm')   gradLm   = opts.gradLm;   else gradLm   = @(rho,m) m./rho; end
if isfield(opts,'funcF') funcF = opts.funcF; else funcF = @(rho) zeros(size(rho)); end
if isfield(opts,'gradF') gradF = opts.gradF; else gradF = @(rho) zeros(size(rho)); end
if isfield(opts,'funcG') funcG = opts.funcG; else funcG = @(rhoend) zeros(size(rhoend)); end
if isfield(opts,'gradG') gradG = opts.gradG; else gradG = @(rhoend) zeros(size(rhoend)); end

% stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 5e3;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end

% parameter for backtracking
if isfield(opts,'stepsize0') stepsize0 = opts.stepsize0; else stepsize0 = 0.1; end
if isfield(opts,'stepmodif') stepmodif = opts.stepmodif; else stepmodif = 0.8; end 
if isfield(opts,'submaxit')  submaxit  = opts.submaxit;  else submaxit  = 5; end

% for solving Poisson equation
q = (2*nt-1:-2:1);
eigenvalueT = ( 2 - 2*cos(pi*q/(2*nt+1)) )/dt/dt;
dA = cell(nt,1);
for i = 1:nt
    dA{i} = decomposition(surf.stiffMatrix + eigenvalueT(i)*surf.massMatrix);
end

dotRho = @(a,b) sum(a.*b.*surf.ptArea,'all')*dt;
dotM = @(a,b) sum(a.*b.*surf.trgArea,'all')*dt;

%% initialization
objArray = zeros(maxit,1);
resArray = zeros(maxit,1);
costArray = cell(maxit,1);
projerrArray = zeros(maxit,1);
stepsizeArray = zeros(maxit,1);

yrho = xrhoNitm;
ym = xmNitm;
% [objNitm,~,~] = compObjGrad(yrho,ym);
stepsize = stepsize0;
wNitm = 1;
%% main iteration
for nit = 1:maxit
    % backtracking
    [objy,gradRho,gradM,~] = compObjGrad(yrho,ym);
    
    subnit = 1;
%     stepsize = stepsize0;
    while subnit < submaxit
        [xrhoNit,xmNit,projerr] = compProj(yrho - stepsize*gradRho, ...
                                           ym - stepsize*gradM);
        xrhoNit = max(xrhoNit,0.01);

        [objNit,~,~,cost] = compObjGrad(xrhoNit,xmNit);
        diffRho = xrhoNit - yrho;
        diffM  = xmNit - ym;
        G = objy + dotRho(gradRho,diffRho) + dotM(gradM,diffM) + ...
            1/(2*stepsize)*(dotRho(diffRho,diffRho)+dotM(diffM,diffM));
        
        if objNit <= G
            break
        end
        subnit = subnit + 1;
        stepsize = stepsize*stepmodif;
    end
    
%     mesh(xrho_k);title(['n=',num2str(nit)]);drawnow;pause(0.005)
    
    % update variables
    wNit = (1 + sqrt(1+4*wNitm^2))/2;
    w = acc*(wNitm-1)/wNit;
%     w = 0;
    drho = xrhoNit-xrhoNitm;
    dm = xmNit -xmNitm;
    yrho = max( xrhoNit + w*drho,0.1);
%     yrho = xrhoNit + w*drho;
    ym  = xmNit  + w*dm;
    
    wNitm = wNit;
    objNitm = objNit;
    xrhoNitm = xrhoNit;
    xmNitm = xmNit;
    stepsize = max(stepsize,1e-5);
    
    res = sqrt(dotRho(drho,drho)+dotM(dm,dm));
    
    objArray(nit) = objNit;
    resArray(nit) = res;
    costArray{nit} = cost;
    projerrArray(nit) = projerr;
    stepsizeArray(nit) = stepsize;
    
    if res < tol
        break
    end
    fprintf('nit%d: %d sub iters, proj err %e, obj %f\n',...
             nit,subnit,projerr,objArray(nit));
    
    if isnan(objNitm) || isinf(objNitm)
        fprintf('blow up at iteration %d with residue %.2e\n',nit,res);
        rho = cat(2,rho0,xrhoNit);
        m = xmNit;
        outs.objArray = objArray(1:nit);
        outs.resArray = resArray(1:nit);
        outs.costArray = costArray(1:nit);
        outs.projerrArray = projerrArray(1:nit);
        outs.stepsizeArray = stepsizeArray(1:nit);
        return
    end
end

%% copy results
rho = cat(2,rho0,xrhoNit);
m = xmNit;
% outs.phi = compPhi(xrhoNit,xmNit);
outs.objArray = objArray(1:nit);
outs.resArray = resArray(1:nit);
outs.costArray = costArray(1:nit);
outs.projerrArray = projerrArray(1:nit);
outs.stepsizeArray = stepsizeArray(1:nit);

%% functions

    function [obj,gradRho,gradM,cost] = compObjGrad(rho,m)
        Frho = funcF(rho(:,1:end-1)).*surf.ptArea;
        Frho(isinf(Frho)|isnan(Frho)) = 0;
        
        Grhoend = funcG(rho(:,end)).*surf.ptArea;
        
        isRhoPtPos = rho > 0;
        rhoend = rho(:,end);
        rhoBar = cat(2,rho0,rho);
        rhoBar = surf.pt2trg*(It(rhoBar));
%         ind = rho > 1e-8;
        isRhoTrgPos = rhoBar > 0;
%         ind = rho ~=0;
        Lrhom = funcL(rhoBar,m).*surf.trgArea;
        
        cost.dyn = dt*(sum(Lrhom(isRhoTrgPos)));
        cost.inter = dt*sum(Frho(:));
        cost.term = sum(Grhoend(:));
        obj = cost.dyn + cost.inter + cost.term;
%         obj = dt*( sum(Lrhom(isRhoTrgPos))+sum(Frho(:)) ) + sum(Grhoend(:));
        
        gradRho = gradLrho(rhoBar,m).*surf.trgArea*dt;
        gradRho(~isRhoTrgPos) = 0;
        gradRho = surf.trg2pt*gradRho;
        gradRho(:,1:end-1) = It(gradRho) + gradF(rho(:,1:end-1)).*surf.ptArea*dt ;
        gradRho(:,end) = gradRho(:,end)/2 + gradG(rhoend).*surf.ptArea;
        gradRho(~isRhoPtPos) = 0;
        
        gradM = gradLm(rhoBar,m).*surf.trgArea*dt;
        gradM = reshape(gradM,[],nDim);
        gradM(~isRhoTrgPos,:) = 0;
        gradM = reshape(gradM,nTrg,nt,nDim);
        
    end

    function [projRho,projM,projerr] = compProj(rho,m)
        %computing the projection
        projRho = rho;
        projM  = m;
        rho = cat(2,rho0,rho);
        
        DtRho = pprDt3(rho);
        rhs = surf.massMatrix*DtRho(:,2:end) + surf.PPPRdivT(m);
        rhs = dctMfg(rhs);
        
        phi = zeros(nPt,nt);
        for j = 1:nt
            phi(:,j) = dA{j}\rhs(:,j);
        end
        phi = idctMfg(phi);
        
        projRho = projRho + pprDt3(phi);
        % projRho = projRho + [pprDt(phi),-phi(:,end)/dt];
%         projRho = max(projRho,0);
        projM  = projM  + surf.PPPRgradT(phi);
        
        %---- proj error
        DtProjRho=pprDt3(cat(2,rho0,projRho));
        projerr = surf.massMatrix*DtProjRho(:,2:end)+surf.divT(projM);
        projerr = max(abs(projerr(:)));
        
    end

    function phi = compPhi(rho,m)
        %computing phi
        [~,Grho,Gm] = compObjGrad(rho,m);
        Grho = pprDt3(Grho);
        % Grho = [-Grho(:,1)/dt, (Grho(:,1:end-1)-Grho(:,2:end))/dt];
        Gm = surf.BxTrg'*squeeze(Gm(:,:,1))...
           + surf.ByTrg'*squeeze(Gm(:,:,2))...
           + surf.BzTrg'*squeeze(Gm(:,:,3));
        b = Grho + Gm;
        
        %%%%%% solve -Delta phi = rho_t + m_x
        B = dctMfg(b);
        phi = zeros(nPt,nt);
        for j = 1:nt
            phi(:,j) = dA{j}\B(:,j);
        end
        phi = idctMfg(phi)*nt;   
    end

    function DtA = Dt(A)
        DtA = (A(:,2:end) - A(:,1:end-1))/dt;
    end

    function ItA = It(A)
        ItA = (A(:,1:end-1) + A(:,2:end))/2;
    end

end