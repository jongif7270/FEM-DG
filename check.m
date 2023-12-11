function [MJ] = check(Ma,N)
[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V);

[c4n,n4e,~,~] = mesh_fem_2d_triangle(0,1,0,1,Ma,Ma,N);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
J=J(1);
I = eye((N+1)*(N+2)/2);
M=I/(V*V');

V=det(V);
M=det(M);
MJ=det(M*J);