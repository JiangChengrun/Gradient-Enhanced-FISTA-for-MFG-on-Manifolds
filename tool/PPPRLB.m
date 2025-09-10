function [Bx, By, Bz] = PPPRLB(mesh)
%% SPRLB Superconvergent Patch Recovery for LB operator
%   RDU = SPRLB(NODE, ELEM, UH) implement Superconvergent Patch Recovery
%   for surface finite element method for Laplace Beltrami Operator
%
%   Reference:
%      Guozhi Dong and Hailong Guo, Polynomial Preserving Recovery
%      of linear finite elements for the Laplace-Beltrami operator on
%      general surfaces, in preparation.
%
%   See also L2PLB, SALB, SPRLB,  WALB.
%
%   Copyright (C) Hailong Guo
%   01/11/2016
%   Reorganized in 04/22/2022

%% Initialization
node = mesh.node;
elem = mesh.elem;
N = size(node, 1);
NT = size(elem, 1);

%% Get Auxiliary structure
T = auxstructure(elem);
elem2edge = double(T.elem2edge);
edge2elem = T.edge2elem;
bdEdge = T.bdEdge;

%% Build point patch
t2v = sparse([1:NT,1:NT,1:NT], double(elem), 1, NT, N);
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
vertn=vertnormal(mesh); % approximated normal vector.

%% Build boundary vertex index
isBdNode = false(N, 1);
bdNode = union(bdEdge(:, 1), bdEdge(:, 2));
isBdNode(bdNode) = true;

%% Least square fitting
% preallocation
Bxij = zeros(10*N, 1);
Byij = zeros(10*N, 1);
Bzij = zeros(10*N, 1);
ii = zeros(10*N, 1);
jj = zeros(10*N, 1);
ind = 0;
for i = 1:N
    nodeStar = find(t2v(:,i));
    if valence(i) < 5
        index = unique(elem2edge(nodeStar,:));
        nodeStar = unique(edge2elem(index,1:2))';
    end
    index = unique(elem(nodeStar, :));
    while isBdNode(i) && length(index) < 7
        index = unique(elem2edge(nodeStar,:));
        nodeStar = unique(edge2elem(index,1:2))';
        index = unique(elem(nodeStar, :));
    end
    pnts = node(index,:);
    NP=size(index,1);
%     uon = surfacedata.unitoutnormal(node(i,:));
    uon=vertn(i,:); % use the approximated normal vector.
    projNode = pnts+ repmat(dot(repmat(node(i,:),NP,1)-pnts, ...
        repmat(uon,NP,1),2), 1, 3).*repmat(uon,NP,1);
    sur_g=dot(pnts-repmat(node(i,:),NP,1), repmat(uon,NP,1),2);% local surface sampled as graph
    if norm(projNode(1,:)-node(i,:)) < eps
        axisu = (projNode(2,:)-node(i,:))./norm(projNode(2,:)-node(i,:));
    else
        axisu = (projNode(1,:)-node(i,:))./norm(projNode(1,:)-node(i,:));
    end
    axisv=cross(uon, axisu);
    pnts = [dot(projNode-repmat(node(i,:),NP,1), repmat(axisu,NP,1), 2), ...
            dot(projNode-repmat(node(i,:),NP,1), repmat(axisv,NP,1), 2)];
    hu = max(abs(pnts(:, 1)));
    hv = max(abs(pnts(:, 2)));
    h = max(hu, hv);
    pnts = pnts./h;
    A = ones(NP, 6);
    A(:, [2, 3]) = pnts;
    A(:, [4, 6]) = pnts.*pnts;
    A(:, 5) = pnts(:, 1).*pnts(:, 2);
    A  =  (A'*A)\A';
    sol_s= A*sur_g;
    rDs =[sol_s(2)/h, sol_s(3)/h];
    J=[eye(2);rDs ];
    J_dash = pinv(J)';
    coord = [axisu', axisv', uon'];
    coe = coord*J_dash;
    Bxij(ind+1:ind+length(index)) = coe(1,1)*A(2, :)/h + ...
        coe(1,2)*A(3,:)/h;
    Byij(ind+1:ind+length(index)) = coe(2,1)*A(2, :)/h + ...
        coe(2,2)*A(3,:)/h;
    Bzij(ind+1:ind+length(index)) = coe(3,1)*A(2, :)/h + ...
        coe(3,2)*A(3,:)/h;
    ii(ind+1:ind+length(index)) = i;
    jj(ind+1:ind+length(index)) = index;
    ind = ind + length(index);
end

%% Construct spare coefficent matrix
Bxij = Bxij(1:ind);
Byij = Byij(1:ind);
Bzij = Bzij(1:ind);
ii = ii(1:ind);
jj = jj(1:ind);
Bx = sparse(ii, jj, Bxij, N, N);
By = sparse(ii, jj, Byij, N, N);
Bz = sparse(ii, jj, Bzij, N, N);
