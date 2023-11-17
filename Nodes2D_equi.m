function [x,y] = Nodes2D_equi(N)
Np = (N+1)*(N+2)/2;
% Create equidistributed nodes on equilateral triangle
L1 = zeros(Np,1); L2 = zeros(Np,1); L3 = zeros(Np,1);
sk = 1;
for n=1:N+1
for m=1:N+2-n
L1(sk) = (n-1)/N; L3(sk) = (m-1)/N;
sk = sk+1;
end
end
L2 = 1.0-L1-L3;
x = -L2+L3; y = (-L2-L3+2*L1)/sqrt(3.0);
return