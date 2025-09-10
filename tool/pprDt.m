function dt_rho = pprDt(rho)
% 对每个空间点独立地, 在时间维做二次拟合恢复导数 ∂t rho
% 输入: rho [Np x (T)], 等间距 t_k = k/(T-1)
% 输出: dt_rho [Np x (T)], 端点用一侧二次, 中间用对称二次
[Np,T] = size(rho);
dt_rho = zeros(Np,T);
if T<3
    % 最少也给个中心差分/一侧差分
    if T==2, dt = 1; dt_rho(:,1) = (rho(:,2)-rho(:,1)); dt_rho(:,2)=dt_rho(:,1); end
    return
end
dt = 1/(T-1);
% 端点：一侧二次 (k=1,2,3 拟合) 与 (k=T-2,T-1,T)
dt_rho(:,1)   = quad_fit_deriv(rho(:,1), rho(:,2), rho(:,3), 0, dt);
dt_rho(:,T)   = quad_fit_deriv(rho(:,T-2), rho(:,T-1), rho(:,T), 2*dt, dt);
% 中间：对称二次 (k-1,k,k+1)
for k=2:T-1
    dt_rho(:,k) = quad_fit_deriv(rho(:,k-1), rho(:,k), rho(:,k+1), dt, dt);
end
end

function d = quad_fit_deriv(y0,y1,y2, t1, dt)
% 拟合 p(t)=a t^2 + b t + c, t∈{0,t1,2dt ?}, 为统一, 取节点等距间隔 dt
% 这里把三个采样点当作 t = 0, dt, 2dt(或 0, t1, 2*dt+t1- dt等价移位) 来做简化
% 我们用标准节点 0,dt,2dt 拟合（平移不影响一阶导数）
Y = [y0, y1, y2]; % Np x 3
% 解 a,b,c: 
% [0^2 0 1; dt^2 dt 1; (2dt)^2 2dt 1] * [a;b;c] = [y0;y1;y2]
M = [0 0 1; dt^2 dt 1; (2*dt)^2 2*dt 1];
coef = (M \ Y.').';   % Np x 3
a = coef(:,1); b = coef(:,2);
d = b;                % p'(0) = b  (中心点用对称窗口时相当于在中点处的导数)
end
