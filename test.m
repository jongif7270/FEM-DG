function err=test(M,N)
xl=0;xr=1;yl=0;yr=1;
[c4n,n4e,ind4e,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,M,M,N);
[ind4e2,~,c4n2] = indexforDG2(xl,xr,yl,yr,M,M,N);

[~,~,~,~,u2,~] = DG3(xl,xr,yl,yr,M,M,N);
[~,u,~,~,~] = FEM(xl,xr,yl,yr,M,M,N);

ue=zeros(size(n4e(:,1),1),(N+1)*(N+2)/2);

[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;
I = eye((N+1)*(N+2)/2);
M=I/(V*V');
err=0;
for j=1:size(n4e,1)
Dx_u = (rx(j)*Dr+sx(j)*Ds)*u(ind4e(j,:));
Dy_u = (ry(j)*Dr+sy(j)*Ds)*u(ind4e(j,:));
Dx_u2 = (rx(j)*Dr+sx(j)*Ds)*u2(ind4e2(j,:));
Dy_u2 = (ry(j)*Dr+sy(j)*Ds)*u2(ind4e2(j,:));
Dex=Dx_u2 - Dx_u;
Dey=Dy_u2 - Dy_u;
err=err+J(j)*(Dex'*M*Dex+Dey'*M*Dey);
end
err=sqrt(err);