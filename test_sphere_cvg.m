% test_sphere_cvg

%function 
close all;
clear;

setpath;
%% 共用参数与代价项
meshName = 'sphere';
[baseSurf.pt,baseSurf.trg] = readOFF(['../Data/',meshName,'.off']);
baseSurf = surfOperators(baseSurf);

% 初末态
makeRho0 = @(pt) exp(-(vecnorm(pt-[cos(pi/3),0,sin(pi/3)],2,2)).^2/0.05) + 0.1;
makeRho1 = @(pt) exp(-(vecnorm(pt-[cos(2*pi/3),0,sin(2*pi/3)],2,2)).^2/0.05) + 0.1;

% 末端 KL 惩罚
buildOpts = @(nt,logrho1) struct( ...
    'plot', 0, ...
    'savegif', 0, ...
    'saveshot', 0, ...
    'nt', nt, ...
    'maxit', 3000, ...
    'tol', 1e-5, ...
    'stepsize0', 8e0, ...
    'stepmodif', 0.5, ...
    'submaxit', 5, ...
    'acc', 0, ...
    'funcL', @(rho,m) sum(m.^2,3)./(2*rho), ...
    'gradLrho', @(rho,m) -sum(m.^2,3)./(2*rho.^2), ...
    'gradLm', @(rho,m) m./rho, ...
    'funcF', @(rho) zeros(size(rho)), ...
    'gradF', @(rho) zeros(size(rho)), ...
    'funcG', @(rhoend) (8e-1)*rhoend.*(log(rhoend)-logrho1), ...
    'gradG', @(rhoend) (8e-1)*(1+log(rhoend)-logrho1) ...
);

%% 时间步长收敛（固定网格）
rho0_base = makeRho0(baseSurf.pt);
rho1_base = makeRho1(baseSurf.pt);
logrho1_base = log(rho1_base);

% 参考解：较细时间步
nt_ref = 256;
opts_ref = buildOpts(nt_ref,logrho1_base);
[rho_ref, ~, ~] = mfgMfFista(baseSurf,rho0_base,opts_ref);

% 被测 nt 列表（从小到大）
nt_list = [16,32,64,128];
time_err = zeros(numel(nt_list),1);
for i = 1:numel(nt_list)
    nt = nt_list(i);
    opts = buildOpts(nt,logrho1_base);
    [rho_nt, ~, ~] = mfgMfFista(baseSurf,rho0_base,opts);
    
    % 与参考解对齐时间：参考解在列索引 1..nt_ref+1，对齐到 nt+1 个时间点
    map_idx = round(linspace(1,nt_ref+1,nt+1));
    rho_ref_sub = rho_ref(:,map_idx);
    
    % L2(时空) 误差：sqrt(∑_{t} ∑_{x} (ρ-ρref)^2 * dA * dt)
    dt = 1/nt;
    diff = rho_nt - rho_ref_sub;
    time_err(i) = sqrt( sum(diff.^2 .* baseSurf.ptArea,'all') * dt );
    fprintf('[Time] nt=%d, L2 spacetime error = %.4e\n',nt,time_err(i));
end

% 时间收敛阶（相邻误差比值的 log2）
time_rate = log(time_err(1:end-1)./time_err(2:end))./log(2);
fprintf('Time-step convergence rates (approx.):\n');
disp(time_rate');

%% 空间网格收敛（固定时间步）
% 生成不同分辨率的球面三角网以做空间收敛
% Fibonacci 采样 + convex hull 得到近似均匀的球面网格
% 若需严格控制网格质量/拓扑，可替换为自有球面网格生成器或外部网格文件。

nt_fixed = 64; % 固定时间步
N_list = [24, 48, 400, 800]; % 顶点数量（越大越细）

% 参考空间解：较细网格 + 同 nt
N_ref = 40000;
fineSurf = makeSphereMesh(N_ref);
fineSurf = surfOperators(fineSurf);
rho0_fine = makeRho0(fineSurf.pt);
rho1_fine = makeRho1(fineSurf.pt);
opts_fine = buildOpts(nt_fixed,log(rho1_fine));
[rho_space_ref, ~, ~] = mfgMfFista(fineSurf,rho0_fine,opts_fine);

space_err = zeros(numel(N_list),1);
for i = 1:numel(N_list)
    Ni = N_list(i);
    surf_i = makeSphereMesh(Ni);
    surf_i = surfOperators(surf_i);
    rho0_i = makeRho0(surf_i.pt);
    rho1_i = makeRho1(surf_i.pt);
    opts_i = buildOpts(nt_fixed,log(rho1_i));
    [rho_i, ~, ~] = mfgMfFista(surf_i,rho0_i,opts_i);
    
    % 将参考解从细网格插值/采样到当前网格（最近邻）
    % 若有 Statistics Toolbox，优先使用 knnsearch；否则退化为 pdist2 最近邻
    ref_idx = nearestOnSphere(surf_i.pt,fineSurf.pt);
    rho_ref_on_i = rho_space_ref(ref_idx,:);
    
    dt = 1/nt_fixed;
    diff = rho_i - rho_ref_on_i;
    space_err(i) = sqrt( sum(diff.^2 .* surf_i.ptArea,'all') * dt );
    fprintf('[Space] N=%d, L2 spacetime error = %.4e\n',Ni,space_err(i));
end

space_rate = log(space_err(1:end-1)./space_err(2:end))./log(2);
fprintf('Space (mesh) convergence rates (approx.):\n');
disp(space_rate');

%% 可视化误差趋势
figure; subplot(1,2,1);
loglog(nt_list,time_err,'o-','LineWidth',1.5); grid on;
xlabel('nt'); ylabel('L2 spacetime error'); title('Time-step convergence');
subplot(1,2,2);
loglog(N_list,space_err,'s-','LineWidth',1.5); grid on;
xlabel('number of vertices'); ylabel('L2 spacetime error'); title('Spatial convergence');

%% 辅助：生成近似均匀球面网格（点+三角）
function surf = makeSphereMesh(N)
    % Fibonacci 采样点
    k = (0:N-1)';
    phi = (1+sqrt(5))/2;
    z = 1 - 2*(k+0.5)/N;
    r = sqrt(max(0,1 - z.^2));
    theta = 2*pi*k/phi;
    x = r.*cos(theta); y = r.*sin(theta);
    P = [x,y,z];
    % 将点外扩一个很小的噪声，避免共球退化
    epsj = 1e-12; P = P + epsj*randn(size(P));
    % 用凸包作为球面三角网
    TR = convhull(P);
    surf.pt = P;
    surf.trg = TR;
end

%% 辅助：最近邻索引（优先 knnsearch）
function idx = nearestOnSphere(Q, P)
    % Q: Mx3 目标点（需在细网格 P 上找最近邻）
    % P: Nx3 参考点
    try
        idx = knnsearch(P,Q);
    catch
        % 退化：用分块 pdist2 避免内存爆
        M = size(Q,1); idx = zeros(M,1);
        blk = 2000;
        for s = 1:blk:M
            e = min(M,s+blk-1);
            D = pdist2(Q(s:e,:),P);
            [~,ii] = min(D,[],2);
            idx(s:e) = ii;
        end
    end
end


