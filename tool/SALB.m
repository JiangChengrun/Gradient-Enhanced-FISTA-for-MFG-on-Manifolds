function [Bx, By, Bz] = SALB(mesh)
%% SALB Simple averaging recovery method for laplace beltrami operator
%
%   RDu = SALB(node, eelm, Du) recover gradient using simple averaging 
%   method for surface finite element method for laplace beltrmai operator
%
%   See also GL2PLB, WALB.
%
%   Copyright (C)  Hailong Guo
%   27/11/2023

%% Initialize
node = mesh.node;
elem = mesh.elem;
NT = size(elem,1);
N = size(node,1);

%% Compute numerical gradient and geometric quantities
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
le1 = sqrt(sum(ve1.^2,2));
le2 = sqrt(sum(ve2.^2,2));
le3 = sqrt(sum(ve3.^2,2));
cp = cross(ve1, ve2, 2);
area =  sqrt(sum(cp.^2,2))/2;
Dlambda(1:NT,:,1) = (repmat(dot(ve1,ve3,2)./le1.^2,1,3).*ve1-ve3).* ...
    repmat(le1.^2./area.^2/4,1,3);
Dlambda(1:NT,:,2) = (repmat(dot(ve2,ve1,2)./le2.^2,1,3).*ve2-ve1).* ...
    repmat(le2.^2./area.^2/4,1,3);
Dlambda(1:NT,:,3) = (repmat(dot(ve2,ve3,2)./le3.^2,1,3).*ve3-ve2).* ...
    repmat(le3.^2./area.^2/4,1,3);
%Du = repmat(uh(elem(:,1)),1,3).*Dlambda(:,:,1) + ...
%    repmat(uh(elem(:,2)),1,3).*Dlambda(:,:,2) + ...
%    repmat(uh(elem(:,3)),1,3).*Dlambda(:,:,3);

%% Recovery gradient
% Note the key part is here for differece between weighted averaging and
% simple averaging. 
area = ones(NT,1);
% dudxArea = area.*Du(:,1);
% dudyArea = area.*Du(:,2);
% dudzArea = area.*Du(:,3);
patchArea = accumarray(elem(:),[area;area;area], [N 1]);
% dudxArea = accumarray(elem(:),[dudxArea;dudxArea;dudxArea],[N 1]);
% dudyArea = accumarray(elem(:),[dudyArea;dudyArea;dudyArea],[N 1]);
% dudzArea = accumarray(elem(:),[dudzArea;dudzArea;dudzArea],[N 1]);
% dudx = dudxArea./patchArea;
% dudy = dudyArea./patchArea;
% dudz = dudzArea./patchArea;
% RDu = [dudx, dudy, dudz];

%% Build point patch
t2v = sparse([1:NT,1:NT,1:NT], elem, 1, NT, N);

%%  Compute the simple average using sparse matrix representation
% preallocation
Bxij = zeros(10*N, 1);
Byij = zeros(10*N, 1);
Bzij = zeros(10*N, 1);
ii = zeros(10*N, 1);
jj = zeros(10*N, 1);
idx = 0;
for i = 1:N
    nodeStar = find(t2v(:,i));
    nV = length(nodeStar)*3;
    Bxij(idx+1:idx+nV) = kron(ones(3,1), area(nodeStar))/patchArea(i).* ...
        reshape(Dlambda(nodeStar, 1, :), nV, 1);
    Byij(idx+1:idx+nV) = kron(ones(3,1), area(nodeStar))/patchArea(i).* ...
        reshape(Dlambda(nodeStar, 2, :), nV, 1);
    Bzij(idx+1:idx+nV) = kron(ones(3,1), area(nodeStar))/patchArea(i).* ...
        reshape(Dlambda(nodeStar, 3, :), nV, 1);
    ii(idx+1:idx+nV) = i;
    jj(idx+1:idx+nV) = reshape(elem(nodeStar, :), nV, 1);
    idx = idx + nV;
end

%% Construct spare coefficent matrix
Bxij = Bxij(1:idx);
Byij = Byij(1:idx);
Bzij = Bzij(1:idx);
ii = ii(1:idx);
jj = jj(1:idx);
Bx = sparse(ii, jj, Bxij, N, N);
By = sparse(ii, jj, Byij, N, N);
Bz = sparse(ii, jj, Bzij, N, N);
