function Dtrho = pprDt3(rho)
%PPRDT3 
[Np,nt] = size(rho);
Dtrho = zeros(Np,nt);
ht = 1/(nt-1);
e = ones(nt,1);

D = spdiags([-e 0*e e],-1:1,nt,nt);
D(1,[1 2 3]) = [-3, 4, -1]; % recovery gradient approximation
D(nt,[nt-2 nt-1 nt]) = [1, -4, 3];% recovery gradient approximation
Bt = D/2/ht;
Dtrho=rho*Bt';