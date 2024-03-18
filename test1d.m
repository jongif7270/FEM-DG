f=@(x) x(:,1).^4+x(:,1).^5;
v=@(x) x(:,1).^3;
% v=@(x) 1;

N=3;
int=2;

[r1D] = JacobiGL(0,0,N);
[int_r1D] = JacobiGL(0,0,int*N);
[V1D] = Vandermonde1D(N,r1D);
[int_V1D]=Vandermonde1D(N,int_r1D);
[iV1D]=Vandermonde1D(int*N,int_r1D);
I1D = eye(size(V1D,1));
iI1D = eye(size(iV1D,1));
A1D=V1D'\int_V1D';
M=I1D/(V1D*V1D');
iM=iI1D/(iV1D*iV1D');

[Dr] = Dmatrix1D(N,r1D,V1D);
[iDr] = Dmatrix1D(int*N,int_r1D,iV1D);

%%integral a) f, b) f' / ans : a) 0.4, b) 2
a=ones(1,N+1)*M*f(r1D);
b=ones(1,N+1)*M*Dr*f(r1D);
%%integral c) vf /int x , d) vf /int o , e) vf /int o cal in / ans : 2/9=0.222
c=ones(1,N+1)*(M.*v(r1D))*f(r1D);
d=ones(1,N+1)*(A1D*(iM.*v(int_r1D))*A1D'*f(r1D));
e=ones(1,N+1)*(A1D*((iM.*v(int_r1D)).*f(int_r1D)')*A1D')*ones(N+1,1);
%%integral g) v'f /int x , h) v'f /int o , i) v'f /int o cal in / ans : 6/7=0.8571
g=ones(1,N+1)*(Dr'*M.*v(r1D))*f(r1D);
h=ones(1,N+1)*(A1D*(iDr'*iM.*v(int_r1D))*A1D'*f(r1D));
i=ones(1,N+1)*(A1D*(iDr'*iM.*v(int_r1D).*f(int_r1D)')*A1D')*ones(N+1,1);

[a b c d e g h i]