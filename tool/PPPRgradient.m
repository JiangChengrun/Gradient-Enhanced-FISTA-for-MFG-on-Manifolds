function [gxN, gyN, gzN] = PPPRgradient(surf,func)
% 构造给PPPRLB用的mesh
mesh.node = surf.pt;     % (nPt x 3)
mesh.elem = surf.trg;    % (nTrg x 3)

% 计算恢复矩阵
[Bx, By, Bz] = PPPRLB(mesh); 
% 向量化一次算完所有时间层
gxN = Bx * func;   % (nPt x nt)
gyN = By * func;   % (nPt x nt)
gzN = Bz * func;   % (nPt x nt)
%% 取某个时间层 j 的节点梯度 (nPt x 3)
% j = round(nt/2);  % 举例取中间时刻
% grad_rho_node = [gxN(:,j), gyN(:,j), gzN(:,j)];

% surf.pt2trg 在 surfOperators.m 里构建 (nTrg x nPt)，做顶点到三角形平均

% grad_rho_trg = cat(3, surf.pt2trg * gxN(:,j), ...
%                        surf.pt2trg * gyN(:,j), ...
%                        surf.pt2trg * gzN(:,j));   % (nTrg x 3)

%nv = surf.normalVec;                          % (nTrg x 3)
% dotn = sum(grad_rho_trg .* nv, 2);            % (nTrg x 1)
% grad_rho_trg_tan = grad_rho_trg - dotn .* nv; % 切向分量

%%%%%%%或改用顶点切向%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vn = vertnormal(mesh);                         % (nPt x 3)
% dotnN = sum(grad_rho_node .* vn, 2);           % (nPt x 1)
% grad_rho_node_tan = grad_rho_node - dotnN .* vn;