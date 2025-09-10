%function test_cvg_parallel
% This function tests the convergence of the MF-Fista algorithm for the
% 2D membrane free-energy functional with PARALLEL COMPUTING
close all;
clear;

setpath;

%% 检查并启动并行计算池
try
    % 检查是否已有并行池
    if isempty(gcp('nocreate'))
        % 启动并行池，使用所有可用核心
        parpool('local');
        fprintf('启动并行计算池，使用 %d 个核心\n', gcp('NumWorkers'));
    else
        fprintf('使用现有并行计算池，%d 个核心\n', gcp('NumWorkers'));
    end
    catch
        warning('无法启动并行计算池，将使用串行计算');
        % 如果无法启动并行池，继续使用串行计算
    end

%% cost funcs 
% klPenalty = 3;
opts.plot = 0;
opts.savegif = 0;
egName = 'ot';
opts.funcL = @(rho,m) sum(m.^2,3)./(rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(rho.^2);
opts.gradLm = @(rho,m) 2*m./rho; 
opts.funcF = @(rho) zeros(size(rho));
opts.gradF = @(rho) zeros(size(rho));

maxits = [497,949,1887,3755];
% opts.maxit = 5000;
opts.tol = 1e-3;
opts.stepsize0 = 1;
opts.stepmodif = 0.5;
opts.submaxit = 5;

%% 并行处理不同网格分辨率
% 定义要测试的网格参数
test_cases = struct();
for i = 1:3
    test_cases(i).nt = 4*2^i;
    test_cases(i).meshsize = 12*2^i;
    test_cases(i).maxit = maxits(i);
    test_cases(i).meshName = ['euc_right_reg',num2str(test_cases(i).meshsize)];
end

fprintf('开始并行计算 %d 个测试用例...\n', length(test_cases));

% 使用 parfor 并行处理
results = cell(length(test_cases), 1);
parfor i = 1:length(test_cases)
    fprintf('Worker %d: 处理网格 %s, nt=%d\n', ...
        labindex, test_cases(i).meshName, test_cases(i).nt);
    
    try
        % 加载网格 - 在 parfor 中需要先声明变量
        meshName = test_cases(i).meshName;
        meshData = load(['../Data/',meshName,'.mat'],'surf');
        surf = meshData.surf;
        
        nPt = size(surf.pt,1);
        nTrg = size(surf.trg,1);
        surf = surfOperators(surf);
        
        % 设置参数
        opts_i = opts;
        opts_i.nt = test_cases(i).nt;
        opts_i.maxit = test_cases(i).maxit;
        dt = 1/opts_i.nt;
        
        % 初始和终止密度
        aroundpt = getaroundpt(surf.pt,surf.trg);
        rho0 = surf.pt(:,1) + 0.5;
        rho1 = ones(nPt,1);
        
        % FISTA 求解
        tic;
        [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts_i);
        solve_time = toc;
        
        % 与真解比较
        W2sq_true = 1/120;
        W2sqerr = abs(output.objArray(end)-W2sq_true);
        
        t = linspace(0,1,opts_i.nt+1);    
        x = surf.pt(:,1);
        sqrtterm = sqrt(2.*x.*t + (0.5*t-1).^2);
        rho_true = (sqrtterm + t-1)./(t.*sqrtterm);
        rho_true(:,1) = rho0; 
        rho_true(:,end) = rho1;
        rhoerr = abs(rho-rho_true);
        rhoerr_l2 = sqrt(sum(rhoerr.^2.*surf.ptArea,'all')*dt);
        rhoerr_l1 = sum(rhoerr.*surf.ptArea,'all')*dt;
        
        t = linspace(dt,1-dt,opts_i.nt);
        x = surf.trgCenter(:,1);
        sqrtterm = sqrt(2.*x.*t + (0.5*t-1).^2);
        flux_true = zeros(nTrg,opts_i.nt,3);
        flux_true(:,:,1) = x./(t.^2) - (3-t)./(2*t.^3).*sqrtterm - (t-1).*(t.^2-4)./(8*t.^3)./sqrtterm - (3*t-4)./(2*t.^3);
        fluxerr = abs(flux-flux_true);
        fluxerr_l2 = sqrt(sum(fluxerr.^2.*surf.trgArea,'all')*dt);
        fluxerr_l1 = sum(fluxerr.*surf.trgArea,'all')*dt;
        
        % 在 parfor 中不能直接调用 save，所以只返回结果
        % 保存结果将在 parfor 循环外进行
        
        % 返回结果结构体
        results{i} = struct(...
            'nt', opts_i.nt, ...
            'meshsize', test_cases(i).meshsize, ...
            'meshName', meshName, ...
            'nPt', nPt, ...
            'nTrg', nTrg, ...
            'solve_time', solve_time, ...
            'W2sqerr', W2sqerr, ...
            'rhoerr_l2', rhoerr_l2, ...
            'rhoerr_l1', rhoerr_l1, ...
            'fluxerr_l2', fluxerr_l2, ...
            'fluxerr_l1', fluxerr_l1, ...
            'output', output, ...
            'rho', rho, ...
            'flux', flux ...
        );
        
        fprintf('Worker %d: 完成 %s, 求解时间 %.2f秒\n', ...
            labindex, meshName, solve_time);
        
    catch
        fprintf('Worker %d: 处理 %s 时出错\n', ...
            labindex, test_cases(i).meshName);
        results{i} = struct('error', '计算失败', 'meshName', test_cases(i).meshName);
    end
end

%% 收集和整理结果
fprintf('\n收集并行计算结果...\n');

% 检查是否有错误并保存结果
valid_results = [];
for i = 1:length(results)
    if ~isfield(results{i}, 'error')
        valid_results = [valid_results; results{i}];
        
        % 在串行循环中保存结果
        meshName = results{i}.meshName;
        nt = results{i}.nt;
        rho = results{i}.rho;
        flux = results{i}.flux;
        output = results{i}.output;
        
        % 保存到文件
        save(['results/',meshName,'_nt',num2str(nt),'_mfp_parallel'], ...
            'rho', 'flux', 'output', 'meshName', 'nt');
        
    else
        fprintf('跳过有错误的结果: %s\n', results{i}.meshName);
    end
end

if isempty(valid_results)
    error('所有测试用例都失败了！');
end

% 按 nt 排序
[~, sort_idx] = sort([valid_results.nt]);
valid_results = valid_results(sort_idx);

% 提取数据用于分析
nt_list = [valid_results.nt]';
nx_list = [valid_results.meshsize]';
W2sqerr_list = [valid_results.W2sqerr]';
rhoerr_l2_list = [valid_results.rhoerr_l2]';
rhoerr_l1_list = [valid_results.rhoerr_l1]';
fluxerr_l2_list = [valid_results.fluxerr_l2]';
fluxerr_l1_list = [valid_results.fluxerr_l1]';
solve_times = [valid_results.solve_time]';

%% 计算收敛率
varerr_l2_list = sqrt(rhoerr_l2_list.^2 + fluxerr_l2_list.^2);
varerr_l1_list = rhoerr_l1_list + fluxerr_l1_list;

fprintf('\n=== 收敛性分析结果 ===\n');
fprintf('时间步数: '); fprintf('%d ', nt_list'); fprintf('\n');
fprintf('网格大小: '); fprintf('%d ', nx_list'); fprintf('\n');
fprintf('求解时间: '); fprintf('%.2fs ', solve_times'); fprintf('\n\n');

if length(nt_list) > 1
    fprintf('W2sqerr:\n'); disp(W2sqerr_list');
    fprintf('收敛阶: '); 
    if length(nt_list) > 1
        rates = log(W2sqerr_list(1:end-1)./W2sqerr_list(2:end))./log(2);
        fprintf('%.2f ', rates); fprintf('\n');
    end
    
    fprintf('\n变量 L2 误差:\n'); disp(varerr_l2_list');
    fprintf('收敛阶: '); 
    if length(nt_list) > 1
        rates = log(varerr_l2_list(1:end-1)./varerr_l2_list(2:end))./log(2);
        fprintf('%.2f ', rates); fprintf('\n');
    end
    
    fprintf('\n变量 L1 误差:\n'); disp(varerr_l1_list');
    fprintf('收敛阶: '); 
    if length(nt_list) > 1
        rates = log(varerr_l1_list(1:end-1)./varerr_l1_list(2:end))./log(2);
        fprintf('%.2f ', rates); fprintf('\n');
    end
end

%% 可视化结果
fig = figure;
set(gcf,'unit','centimeters','position',[10 5 15 6]);

% 目标函数误差
subplot(1,2,1);
set(gca,'Position',[0.1,0.2,0.35,0.7]);
loglog(nt_list, W2sqerr_list, 'ro-', 'linewidth', 1.5, 'markersize', 7);
grid on;
xlabel('nt', 'fontsize', 10, 'fontweight', 'bold');
ylabel('objective error', 'fontsize', 10, 'fontweight', 'bold');
title('目标函数误差收敛性');

% 变量 L2 误差
subplot(1,2,2);
set(gca,'Position',[0.6,0.2,0.35,0.7]);
loglog(nt_list, varerr_l2_list, 'bs--', 'linewidth', 1.5, 'markersize', 7);
grid on;
xlabel('nt', 'fontsize', 10, 'fontweight', 'bold');
ylabel('L2 norm of variable error', 'fontsize', 10, 'fontweight', 'bold');
title('变量 L2 误差收敛性');

% 保存图片
exportgraphics(fig, 'results/cvg_parallel.eps');

% 创建额外的分析图
figure;
subplot(2,2,1);
loglog(nt_list, rhoerr_l2_list, 'b-o', 'linewidth', 1.5);
grid on; xlabel('nt'); ylabel('rho L2 error'); title('密度 L2 误差');

subplot(2,2,2);
loglog(nt_list, fluxerr_l2_list, 'r-s', 'linewidth', 1.5);
grid on; xlabel('nt'); ylabel('flux L2 error'); title('通量 L2 误差');

subplot(2,2,3);
semilogx(nt_list, solve_times, 'g-^', 'linewidth', 1.5);
grid on; xlabel('nt'); ylabel('solve time (s)'); title('求解时间');

subplot(2,2,4);
loglog(nx_list, varerr_l2_list, 'm-d', 'linewidth', 1.5);
grid on; xlabel('mesh size'); ylabel('variable L2 error'); title('空间收敛性');

%% 保存结果
save(['results/euc_cvg_parallel_', num2str(max(maxits))], ...
    'valid_results', 'nt_list', 'nx_list', 'W2sqerr_list', ...
    'rhoerr_l2_list', 'rhoerr_l1_list', 'fluxerr_l2_list', ...
    'fluxerr_l1_list', 'solve_times', 'varerr_l2_list', 'varerr_l1_list');

fprintf('\n结果已保存到 results/ 目录\n');
fprintf('并行计算完成！\n');

%% 清理并行池（可选）
% 如果不再需要并行计算，可以关闭
% if ~isempty(gcp('nocreate'))
%     delete(gcp('nocreate'));
%     fprintf('已关闭并行计算池\n');
% end
