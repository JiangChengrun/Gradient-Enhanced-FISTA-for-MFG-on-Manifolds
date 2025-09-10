function Dtrho = pprDt2(rho)
%PPRDT2 此处显示有关此函数的摘要
%   此处显示详细说明
[Np,T] = size(rho);
Dtrho = zeros(Np,T);
dt = 1/(T-1);
Dtrho(:,2:T-1) = (rho(:,3:T) - rho(:,1:T-2)) / (2*dt);                 % 中心二阶
Dtrho(:,1)     = (-3*rho(:,1) + 4*rho(:,2) - rho(:,3)) / (2*dt);       % 前向二阶
Dtrho(:,T)     = ( 3*rho(:,T) - 4*rho(:,T-1) + rho(:,T-2)) / (2*dt);   % 后向二阶

