
% test_sphere_h_cvg.m
% 网格尺寸收敛率（h-收敛）测试：在球面上用 MF-FISTA 求解曲面 MFG，
% 逐次细化三角网格，固定时间离散 nt，比较不同网格下的解与“最细网格参考解”之间的误差，
% 并在 loglog 图上给出参考斜率线以直观看出收敛阶。
%
% 依赖：本仓库中的 SurfOp/ 与 FISTA/ 代码（surfOperators, readOFF, mfgMfFista 等）。
% 建议将本文件放在仓库根目录下（与 test_sphere.m 同级）。
%
% 作者：ChatGPT（根据你的项目结构自动适配）
% --------------------------------------------------------------
close all; clear; clc;
setpath;  % 将当前仓库加入 MATLAB path（见 setpath.m）

%% ----------------- 固定的时间离散与 FISTA 参数 -----------------
opts = struct();
opts.plot = 0;
opts.savegif = 0;
opts.saveshot = 0;

opts.nt = 32;                % 固定时间步数（可按需调整为 16/32/64）
opts.maxit = 2500;           % FISTA 最大迭代
opts.tol = 1e-2;             % 停止阈值
opts.stepsize0 = 8e0;
opts.stepmodif = 0.5;
opts.submaxit = 5;
opts.acc = 0;

% 动力学项（与 test_cvg / test_sphere 一致的最常用形式）
opts.funcL    = @(rho,m) sum(m.^2,3)./(rho);                   % \int |m|^2 / rho
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(rho.^2);
opts.gradLm   = @(rho,m) 2*m./rho;

% 交互/终端项：这里设为 0，便于聚焦网格误差（如需势能，可改写）
opts.funcF = @(rho) zeros(size(rho));
opts.gradF = @(rho) zeros(size(rho));
opts.funcG = @(rho) zeros(size(rho));
opts.gradG = @(rho) zeros(size(rho));

%% ----------------- 初、末态设定（可与 test_sphere 一致） -----------------
% 使用两个球面高斯作为测试（可按需更改）
makeRho0 = @(pt) exp(-(vecnorm(pt-[cos(pi/3),0,sin(pi/3)],2,2)).^2/0.05) + 0.1;
makeRho1 = @(pt) exp(-(vecnorm(pt-[cos(pi/3),0,-sin(pi/3)],2,2)).^2/0.05) + 0.1;

%% ----------------- 读取初始网格（单位球） -----------------
[baseSurf.pt, baseSurf.trg] = readOFF(['../Data/','sphere','.off']);
baseSurf = surfOperators(baseSurf);

%% ----------------- 网格细化层数与参考层 -----------------
nLevels = 3;            % 细化层数（包含最粗层）。建议 3~4 层；过多层计算会较重。
refLevel = nLevels;     % 以“最细层”为参考解层

% 收敛曲线数据初始化
h_list          = zeros(nLevels,1);   % 代表性网格尺寸（sqrt(mean(trgArea))）
ndof_list       = zeros(nLevels,1);   % 顶点数（自由度）
obj_list        = zeros(nLevels,1);   % 目标函数值 J_h
rho_err_list    = nan(nLevels,1);     % L2(时空)-误差 ||rho_h - rho_ref||
flux_err_list   = nan(nLevels,1);     % L2(时空)-误差 ||m_h   - m_ref||

% 存储每一层的解（用于后续与参考层比较）
RHO_all  = cell(nLevels,1);   % rho: (nPt x (nt+1))
M_all    = cell(nLevels,1);   % m:   (nTrg x nt x 3)
SURF_all = cell(nLevels,1);   % 对应的几何/算子

% 逐层细化并求解
surf = baseSurf;
for lev = 1:nLevels
    if lev > 1
        surf = refineSphereMeshOnce(surf);   % 边中点细化并投影回单位球
        surf = surfOperators(surf);
    end
    SURF_all{lev} = surf;
    nPt  = size(surf.pt,1);
    nTrg = size(surf.trg,1);
    ndof_list(lev) = nPt;

    % 网格尺寸 h 的一个经验取法：sqrt(mean(trgArea))
    h_list(lev) = sqrt(mean(surf.trgArea));

    % 初/末态在该层网格上的取值
    rho0 = makeRho0(surf.pt);
    rho1 = makeRho1(surf.pt);

    % ---------- FISTA 求解 MFG ----------
    % 你可以自由切换 mfgMfFista（只给 rho0）或 mfpMfFista（同时给 rho0,rho1）
    % 下面默认用 mfpMfFista，有末态约束；若想与 test_sphere 完全一致也可换成 mfgMfFista。
    fprintf('Level %d/%d: nPt=%d, nTrg=%d ...\n', lev, nLevels, nPt, nTrg);
    % [rho, m, outs] = mfpMfFista(surf, rho0, rho1, opts);
    [rho, m, outs] = GEmfgMfFista(surf, rho0, opts);

    % 记录
    RHO_all{lev} = rho;
    M_all{lev}   = m;

    % 动态项 + 交互 + 终端的总目标值（与 test_sphere 用法一致）
    dt   = 1/opts.nt;
    rhos = surf.pt2trg*(rho(:,1:end-1)+rho(:,2:end))/2;    % 三角片上时间中心的 rho
    Jdyn = sum(surf.trgArea .* opts.funcL(rhos, m), 'all') * dt;
    Jint = sum(surf.ptArea .* opts.funcF(rho(:,2:end-1)), 'all') * dt;
    Jter = sum(surf.ptArea .* opts.funcG(rho(:,end)), 'all');
    obj_list(lev) = Jdyn + Jint + Jter;   % 也可用 outs.objArray(end)

end

%% ----------------- 与最细层参考解比较（零阶近邻延拓） -----------------
% 为避免实现复杂的跨网格插值，这里使用“最近邻点/片元”将粗网格解延拓到最细网格，
% 然后在最细网格上做 L2(时空) 误差；对 m 使用三角片心的最近邻。
surfRef = SURF_all{refLevel};
rhoRef  = RHO_all{refLevel};    % (nPt_ref x (nt+1))
mRef    = M_all{refLevel};      % (nTrg_ref x nt x 3)
dt      = 1/opts.nt;

% 预先计算三角片心
bcRef = ( surfRef.pt(surfRef.trg(:,1),:) + ...
          surfRef.pt(surfRef.trg(:,2),:) + ...
          surfRef.pt(surfRef.trg(:,3),:) )/3;

for lev = 1:nLevels
    if lev == refLevel
        rho_err_list(lev)  = 0;
        flux_err_list(lev) = 0;
        continue;
    end
    surfC = SURF_all{lev};
    rhoC  = RHO_all{lev};    % (nPt_c x (nt+1))
    mC    = M_all{lev};      % (nTrg_c x nt x 3)

    % 1) 将粗层顶点的 rho 延拓到“参考层”的顶点：最近邻点
    idxPt = nearestNeighborIndex(surfC.pt, surfRef.pt);   % size: nPt_ref x 1
    rhoC_on_ref = rhoC(idxPt, :);                          % (nPt_ref x (nt+1))

    % 2) 将粗层三角片的 m 延拓到“参考层”的三角片：用片心最近邻
    bcC = ( surfC.pt(surfC.trg(:,1),:) + ...
            surfC.pt(surfC.trg(:,2),:) + ...
            surfC.pt(surfC.trg(:,3),:) )/3;
    idxTrg = nearestNeighborIndex(bcC, bcRef);            % size: nTrg_ref x 1
    % 重新索引到参考层的三角数量
    mC_on_ref = mC(idxTrg, :, :);                          % (nTrg_ref x nt x 3)

    % 3) 计算 L2(时空) 误差
    rhoDiff = rhoC_on_ref - rhoRef;                        % (nPt_ref x (nt+1))
    rho_err_list(lev) = sqrt( sum( (rhoDiff.^2) .* surfRef.ptArea, 'all') * dt );

    mDiff = mC_on_ref - mRef;                              % (nTrg_ref x nt x 3)
    mDiff2 = sum(mDiff.^2, 3);                             % |m|^2 (nTrg_ref x nt)
    flux_err_list(lev) = sqrt( sum( mDiff2 .* surfRef.trgArea, 'all') * dt );
end

%% ----------------- 目标函数的收敛（相对最细层） -----------------
obj_ref = obj_list(refLevel);
obj_err = abs(obj_list - obj_ref);

%% ----------------- 拟合斜率并画 loglog 收敛图 -----------------
% 常规习惯是用“从粗到细的 h”；确保 h_list 单调递减
[h_list, perm] = sort(h_list(:), 'descend');
rho_err_plot   = rho_err_list(perm);
flux_err_plot  = flux_err_list(perm);
obj_err_plot   = obj_err(perm);
ndof_plot      = ndof_list(perm);

% 去掉参考层的 0 误差点以作线性拟合（避免 -Inf）
mask = rho_err_plot > 0;
p_rho  = polyfit(log(h_list(mask)),  log(rho_err_plot(mask)), 1);  slope_rho  = p_rho(1);
mask2 = flux_err_plot > 0;
p_flux = polyfit(log(h_list(mask2)), log(flux_err_plot(mask2)),1); slope_flux = p_flux(1);
mask3 = obj_err_plot > 0;
p_obj  = polyfit(log(h_list(mask3)), log(obj_err_plot(mask3)), 1); slope_obj  = p_obj(1);

% 参考线通过最细层的点： y_ref = C * h^{slope}
C_rho  = rho_err_plot(find(mask,1,'last')) / h_list(find(mask,1,'last'))^slope_rho;
C_flux = flux_err_plot(find(mask2,1,'last'))/ h_list(find(mask2,1,'last'))^slope_flux;
C_obj  = obj_err_plot(find(mask3,1,'last')) / h_list(find(mask3,1,'last'))^slope_obj;

% 图 1：rho 的 L2(时空) 误差
figure; loglog(h_list, rho_err_plot, 'o-','linewidth',1.6,'markersize',6); hold on;
loglog(h_list, C_rho*h_list.^slope_rho, 'k--','linewidth',1.2);
grid on; xlabel('h'); ylabel('||\rho_h - \rho_{ref}||_{L^2(S\times[0,1])}');
title(sprintf('rho: 拟合斜率 \\approx %.2f', slope_rho));
legend({'误差','参考斜率线'}, 'location','southwest');
set(gcf,'color','w');

if ~exist('results','dir'); mkdir results; end
exportgraphics(gcf, fullfile('results','sphere_h_cvg_rho.png'));


% 图 2：m 的 L2(时空) 误差
figure; loglog(h_list, flux_err_plot, 's-','linewidth',1.6,'markersize',6); hold on;
loglog(h_list, C_flux*h_list.^slope_flux, 'k--','linewidth',1.2);
grid on; xlabel('h'); ylabel('||m_h - m_{ref}||_{L^2(S\times[0,1])}');
title(sprintf('flux: 拟合斜率 \\approx %.2f', slope_flux));
legend({'误差','参考斜率线'}, 'location','southwest');
set(gcf,'color','w');

if ~exist('results','dir'); mkdir results; end
exportgraphics(gcf, fullfile('results','sphere_h_cvg_flux.png'));


% 图 3：目标函数 J 的收敛
figure; loglog(h_list, obj_err_plot, 'd-','linewidth',1.6,'markersize',6); hold on;
loglog(h_list, C_obj*h_list.^slope_obj, 'k--','linewidth',1.2);
grid on; xlabel('h'); ylabel('|J_h - J_{ref}|');
title(sprintf('目标函数: 拟合斜率 \\approx %.2f', slope_obj));
legend({'误差','参考斜率线'}, 'location','southwest');
set(gcf,'color','w');

if ~exist('results','dir'); mkdir results; end
exportgraphics(gcf, fullfile('results','sphere_h_cvg_obj.png'));


%% ----------------- 实用函数：一次细化（边中点 + 单位球投影） -----------------
function surf2 = refineSphereMeshOnce(surf1)
    % 输入：surf1.pt (n x 3), surf1.trg (m x 3)
    % 输出：surf2 的 pt,trg（每个旧三角细分为 4 个新三角）
    pt  = surf1.pt;
    trg = surf1.trg;
    nPt = size(pt,1);
    nTrg= size(trg,1);

    % 用 map 记录已经生成的中点，键为 'i_j' (i<j)
    edgeMap = containers.Map('KeyType','char','ValueType','int32');
    newPt = pt;  % 先复制旧点
    % 新三角集合（预分配近似长度）
    newTrgs = zeros(nTrg*4, 3);
    tcount = 0;

    for T = 1:nTrg
        i = trg(T,1); j = trg(T,2); k = trg(T,3);
        % 三条边的中点索引
        ij = edgeMidIndex(i,j);
        jk = edgeMidIndex(j,k);
        ki = edgeMidIndex(k,i);
        % 四个子三角
        tcount = tcount + 1; newTrgs(tcount,:) = [i,  ij, ki];
        tcount = tcount + 1; newTrgs(tcount,:) = [ij, j,  jk];
        tcount = tcount + 1; newTrgs(tcount,:) = [ki, jk, k ];
        tcount = tcount + 1; newTrgs(tcount,:) = [ij, jk, ki];
    end
    newTrgs = newTrgs(1:tcount,:);
    surf2.pt  = newPt;
    surf2.trg = newTrgs;

    function idx = edgeMidIndex(a,b)
        if a>b, t=a; a=b; b=t; end
        key = sprintf('%d_%d',a,b);
        if isKey(edgeMap,key)
            idx = edgeMap(key);
            return;
        end
        % 生成新点：边中点并重新投影到单位球
        mid = (pt(a,:) + pt(b,:))/2;
        mid = mid / norm(mid,2);
        newPt = [newPt; mid];
        idx = size(newPt,1);
        edgeMap(key) = idx;
    end
end

%% ----------------- 实用函数：最近邻索引（无工具箱依赖） -----------------
function idx = nearestNeighborIndex(srcPts, tgtPts)
    % 对每个 tgtPts(q,:) 找到 srcPts 中最近的点索引
    % 优先使用 knnsearch（若可用），否则退化为分块暴力搜索。
    try
        idx = knnsearch(srcPts, tgtPts);
        return;
    catch
        % 分块计算，避免一次性构造巨大距离矩阵
        nSrc = size(srcPts,1);
        nTgt = size(tgtPts,1);
        blk  = 5000; % 每块处理这么多目标点，按需调大/调小
        idx = zeros(nTgt,1);
        for s = 1:blk:nTgt
            e = min(s+blk-1, nTgt);
            B = tgtPts(s:e,:);                  % (b x d)
            % 计算 (b x nSrc) 的距离平方
            % dist^2 = |a|^2 + |b|^2 - 2 a·b
            a2 = sum(srcPts.^2,2)';             % (1 x nSrc)
            b2 = sum(B.^2,2);                   % (b x 1)
            D2 = bsxfun(@plus, a2, b2) - 2*(B*srcPts');  %#ok<BSXFUN>
            [~, ii] = min(D2, [], 2);
            idx(s:e) = ii;
        end
    end
end
