function mesh = preprocessTS(mesh)
%% PREPROCESS Preprocess for gradient and divergence
%
%   Copyright (C) Hailong Guo
%   23/09/2022

%% Initialize
nt = mesh.nt;
%hx = mesh.hx;
ht = 1/(nt-1);
e = ones(nt,1);

%% Construct the recovery matrix in time and space
if strcmp(mesh.recoverymethod, 'PPPR')
    %Differentiation matrix in time
    D = spdiags([-e 0*e e],-1:1,nt,nt);
    D(1,[1 2 3]) = [-3, 4, -1]; % recovery gradient approximation
    D(nt,[nt-2 nt-1 nt]) = [1, -4, 3];% recovery gradient approximation
    Bt = D/2/ht;
    % Differentiation matrix in space
    [Bx, By, Bz] = PPPRLB(mesh);
elseif strcmp(mesh.recoverymethod, 'SA')
    %Differentiation matrix in time
    Bt = 0.5*spdiags([-e 0*e e],-1:1,nt,nt)/ht;
    Bt(1, 1) = -1/ht;
    Bt(1, 2) = 1/ht;
    Bt(nt, nt-1) = -1/ht;
    Bt(nt, nt) = 1/ht;
    % Differentiation matrix in space
    [Bx, By, Bz] = SALB(mesh);
else
    error('The recovery method is not implemented')
end

%% Parse the output
mesh.Bx = Bx;
mesh.By = By;
mesh.Bz = Bz;
mesh.Bt = Bt;
