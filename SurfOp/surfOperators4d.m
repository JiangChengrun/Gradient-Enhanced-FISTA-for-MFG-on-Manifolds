function surf = surfOperators4d(surf)
    % This function computes the surface operators for a given surface
    % 1. Compute the normal vector, area, and center of each triangle
    % 2. Compute the gradient and divergence operators
    % 3. Compute the mass and stiffness matrices
    % 4. Compute the gradient and divergence operators in tangential direction
    % 5. Compute the inverse distance weighting operators
    % 6. Compute the second fundamental form
    

if nargin<1
    disp('Parameter:  Surf = SurfOperators(Surf)');
    return;
end

pt = surf.pt;
trg = surf.trg;

nPt = length(pt);
nTrg = length(trg);
trgCenter = zeros(nTrg,4);
trgArea = zeros(nTrg,1);
ptArea = zeros(nPt,1);
grad1 = sparse(nTrg,nPt);
grad2 = sparse(nTrg,nPt);
grad3 = sparse(nTrg,nPt);
grad4 = sparse(nTrg,nPt);

for i=1:size(trg,1)
    p1 = trg(i,1);  p2 = trg(i,2); p3 = trg(i,3);
    v1 = pt(p1,:);  v2 = pt(p2,:); v3 = pt(p3,:);
    e1 = v3-v2;     e2 = v1-v3;    e3 = v2-v1;
    l1 = norm(e1);  l2 = norm(e2); l3 = norm(e3);
    
    s  =  (l1+l2+l3)*0.5;
    trgArea(i) = sqrt(s*(s-l1)*(s-l2)*(s-l3));    
    trgCenter(i,:) = mean([v1;v2;v3]);
%     v12 = v2-v1; v31 = v1-v3;
%     trgCenter(i,:) = mean([v1;v2;v3]);
%     trgArea(i) = 1/2*norm(n);
    
    ptArea(p1) = ptArea(p1) + trgArea(i)/3;
    ptArea(p2) = ptArea(p2) + trgArea(i)/3;
    ptArea(p3) = ptArea(p3) + trgArea(i)/3;
    
    
    shapegrad1 = [-1,-1]*[dot(e3,e3), dot(e3,-e2);
                          dot(e3,-e2),dot(e2,e2) ]^(-1)...
                        *[e3;-e2];
    shapegrad2 = [-1,-1]*[dot(e1,e1), dot(e1,-e3);
                          dot(e1,-e3),dot(e3,e3) ]^(-1)...
                        *[e1;-e3];
    shapegrad3 = [-1,-1]*[dot(e2,e2), dot(e2,-e1);
                          dot(e2,-e1),dot(e1,e1) ]^(-1)...
                        *[e2;-e1];
    grad1(i,[p1 p2 p3]) = [shapegrad1(1),shapegrad2(1),shapegrad3(1)];
    grad2(i,[p1 p2 p3]) = [shapegrad1(2),shapegrad2(2),shapegrad3(2)];
    grad3(i,[p1 p2 p3]) = [shapegrad1(3),shapegrad2(3),shapegrad3(3)];
    grad4(i,[p1 p2 p3]) = [shapegrad1(4),shapegrad2(4),shapegrad3(4)];
%     shapegrad = [dot(v12,v12),dot(v12,-v31);
%                  dot(v12,-v31),dot(v31,v31)]^(-1) * [v12;-v31];
%     grad1(i,[p1 p2 p3]) = [-shapegrad(1,1)- shapegrad(2,1), ...
%                             shapegrad(1,1), shapegrad(2,1)];
%     grad2(i,[p1 p2 p3]) = [-shapegrad(1,2)- shapegrad(2,2), ...
%                             shapegrad(1,2), shapegrad(2,2)];
%     grad3(i,[p1 p2 p3]) = [-shapegrad(1,3)- shapegrad(2,3), ...
%                             shapegrad(1,3), shapegrad(2,3)];
%     grad4(i,[p1 p2 p3]) = [-shapegrad(1,4)- shapegrad(2,4), ...
%                             shapegrad(1,4), shapegrad(2,4)];
    
end
% surf.normalVec = normalVec;
surf.ptArea = ptArea;
surf.trgArea = trgArea;
surf.trgCenter = trgCenter;
surf.grad1 = grad1;
surf.grad2 = grad2;
surf.grad3 = grad3;
surf.grad4 = grad4;

%---- trg2pt pt2trg 
I = repmat((1:nTrg)',1,3);
J = surf.trg;
surf.pt2trg = sparse(I(:),J(:),1/3*ones(3*nTrg,1),nTrg,nPt);
surf.trg2pt = surf.pt2trg';

%---- massmat,stiffmat
surf.massMatrix = spdiags(ptArea,0,nPt,nPt);
surf.stiffMatrix = grad1'*(trgArea.*grad1)...
                  +grad2'*(trgArea.*grad2)...
                  +grad3'*(trgArea.*grad3)...
                  +grad4'*(trgArea.*grad4);



%---- grad, div
surf.grad = @(u) [grad1*u,grad2*u,grad3*u,grad4*u];
surf.div = @(V) -( grad1'*(trgArea.*V(:,1)) ...
                  +grad2'*(trgArea.*V(:,2)) ...
                  +grad3'*(trgArea.*V(:,3)) ...
                  +grad4'*(trgArea.*V(:,4)) );
surf.gradT = @(u) full(cat(3,(grad1*u),(grad2*u),(grad3*u),(grad4*u)));
surf.divT = @(V) -( grad1'*(trgArea.*V(:,:,1)) ...
                   +grad2'*(trgArea.*V(:,:,2)) ...
                   +grad3'*(trgArea.*V(:,:,3)) ...
                   +grad4'*(trgArea.*V(:,:,4)) );

end
