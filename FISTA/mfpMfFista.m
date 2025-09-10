function [rho,m,outs] = mfpMfFista(surf,rho0,rho1,opts)
    % This function computes the MF-Fista algorithm for minimizing the
    % 2D membrane free-energy functional
    % Input: surf = surface operator
    %        rho0 = initial value of rho
    %        rho1 = final value of rho
    %        opts = options
    % Output: rho = final value of rho
    %         m = final value of m
    %         outs = output structure       
%% parameters
nPt = size(surf.pt,1);
nTrg = size(surf.trg,1);
nDim = size(surf.pt,2);

if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
ntm = nt-1;
dt = 1/nt;

% inital value of rho and m
if isfield(opts,'rho') xrhoNitm = opts.rho; else xrhoNitm = ones(nPt,ntm); end
if isfield(opts,'m')   xmNitm  = opts.m;    else xmNitm  = zeros(nTrg,nt,nDim); end

% functions and gradients
if isfield(opts,'funcL') funcL = opts.funcL; else funcL = @(rho,m) sum(m.^2,3)./rho; end
if isfield(opts,'gradLrho') gradLrho = opts.gradLrho; else gradLrho = @(rho,m) -sum(m.^2,3)./(rho.^2); end
if isfield(opts,'gradLm')   gradLm   = opts.gradLm;   else gradLm   = @(rho,m) 2*m./rho; end
if isfield(opts,'funcF') funcF = opts.funcF; else funcF = @(rho) zeros(size(rho)); end
if isfield(opts,'gradF') gradF = opts.gradF; else gradF = @(rho) zeros(size(rho)); end

% stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 5e3;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end

% parameter for backtracking
if isfield(opts,'stepsize0') stepsize0 = opts.stepsize0; else stepsize0 = 0.1; end
if isfield(opts,'stepmodif') stepmodif = opts.stepmodif; else stepmodif = 0.8; end 
if isfield(opts,'submaxit')  submaxit  = opts.submaxit;  else submaxit  = 5; end

% for solving Poisson equation
eigenvalueT = ( 2 - 2*cos(pi*(0:nt-1)/nt) )/dt/dt;
dA = cell(nt,1);
for i = 1:nt
    dA{i} = decomposition(surf.stiffMatrix + eigenvalueT(i)*surf.massMatrix);
end

dotRho = @(a,b) sum(a.*b.*surf.ptArea,'all')*dt;
dotM = @(a,b) sum(a.*b.*surf.trgArea,'all')*dt;

%% initialization
objArray = zeros(maxit,1);
resArray = zeros(maxit,1);
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
    [objy,gradRho,gradM] = compObjGrad(yrho,ym);
    
    subnit = 1;
%     stepsize = stepsize0;
    while subnit < submaxit
        [xrhoNit,xmNit,projerr] = compProj(yrho - stepsize*gradRho, ym - stepsize*gradM);
        
        [objNit,~,~] = compObjGrad(xrhoNit,xmNit);
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
    w = (wNitm-1)/wNit;
%     w = 0;
    drho = xrhoNit-xrhoNitm;
    dm = xmNit -xmNitm;
    yrho = max( xrhoNit + w*drho,0.1);
%     yrho = xrho_k + w_k*drho;
    ym  = xmNit  + w*dm;
    
    wNitm = wNit;
    objNitm = objNit;
    xrhoNitm = xrhoNit;
    xmNitm = xmNit;
    stepsize = max(stepsize,1e-5);
    
    res = sqrt(dotRho(drho,drho)+dotM(dm,dm));
    
    objArray(nit) = objNit;
    resArray(nit) = res;
    projerrArray(nit) = projerr;
    stepsizeArray(nit) = stepsize;
    
%     if res < tol
%         break
%     end
    fprintf('nit%d: %d sub iters, proj err %e, obj %f\n',...
             nit,subnit,projerr,objArray(nit));
    
    if isnan(objNitm)
        fprintf('blow up at iteration %d with residue %.2e\n',nit,res);
        rho = cat(2,rho0,xrhoNit,rho1);
        m = xmNit;
        outs.objArray = objArray(1:nit);
        outs.resArray = resArray(1:nit);
        outs.projerrArray = projerrArray(1:nit);
        outs.stepsizeArray = stepsizeArray(1:nit);
        return
    end
end

%% copy results
rho = cat(2,rho0,xrhoNit,rho1);
m = xmNit;
outs.objArray = objArray(1:nit);
outs.resArray = resArray(1:nit);
outs.projerrArray = projerrArray(1:nit);
outs.stepsizeArray = stepsizeArray(1:nit);

%% functions

    function [obj,gradRho,gradM] = compObjGrad(rho,m)
        Frho = funcF(rho).*surf.ptArea;
        Frho(isinf(Frho)|isnan(Frho)) = 0;
        
        rhoBar = cat(2,rho0,rho,rho1);
        rhoBar = surf.pt2trg*(It(rhoBar));
%         ind = rho > 1e-8;
        ind = rhoBar > 0;
%         ind = rho ~=0;        
        Lrhom = funcL(rhoBar,m).*surf.trgArea;
        
        obj =  (sum(Frho(:))+sum(Lrhom(ind)))*dt;
        
        gradRho = gradLrho(rhoBar,m).*surf.trgArea;
        gradRho(~ind) = 0;
        gradRho = It(surf.trg2pt*gradRho) + gradF(rho).*surf.ptArea;
        gradM = gradLm(rhoBar,m).*surf.trgArea;

    end

    function [projRho,projM,projerr] = compProj(rho,m)
        projRho = rho;
        projM  = m;
        rho = cat(2,rho0,rho,rho1);
        
        rhs = surf.massMatrix*Dt(rho) + surf.divT(m);
        rhs = dct(rhs,[],2);
        
        phi = zeros(nPt,nt);
        for j = 1:nt
            phi(:,j) = dA{j}\rhs(:,j);
        end
        phi = idct(phi,[],2);
%         %---- linear solver error
%         phi_tt = 2*phi;
%         phi_tt(:,[1 end]) = phi(:,[1 end]);
%         phi_tt(:,1:end-1) = phi_tt(:,1:end-1) - phi(:,2:end);
%         phi_tt(:,2:end) = phi_tt(:,2:end) - phi(:,1:end-1);
%         phi_tt = surf.ptArea.*phi_tt/dt/dt;
%         lhs = phi_tt - surf.stiffMatrix*phi;
%         rhs = surf.ptArea.*Dt(rho) + surf.divT(m);
%         linsolerr = lhs - rhs;
%         linsolerr = max(abs(linsolerr(:)));
%         %----
        
        projRho = projRho + Dt(phi);
        projM  = projM  + surf.gradT(phi);
        
        %---- proj error
        projerr = surf.massMatrix*Dt(cat(2,rho0,projRho,rho1))+surf.divT(projM);
        projerr = max(abs(projerr(:)));
        
    end

    function DtA = Dt(A)
        DtA = (A(:,2:end) - A(:,1:end-1))/dt;
    end

    function ItA = It(A)
        ItA = (A(:,1:end-1) + A(:,2:end))/2;
    end

end