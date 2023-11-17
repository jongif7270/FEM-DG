function [A,V,Dr,Ds,u,ux,uy] = DG(xl,xr,yl,yr,Mx,My,N)
[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,inddb,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);
[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V);
b = zeros(size(ind4e(:),1),1);
u = b;
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;
ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));

[B]=lineint2(xl,xr,yl,yr,Mx,My,N);

I = eye((N+1)*(N+2)/2);
xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

Kr=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
M=J(j)*I/(V*V');
Srr=J(j)*(V\Dr)'*(V\Dr);
Srs=J(j)*(V\Dr)'*(V\Ds);
Ssr=J(j)*(V\Ds)'*(V\Dr);
Sss=J(j)*(V\Ds)'*(V\Ds);

K=J(j)*((rx(j)^2+ry(j)^2)*Srr+(rx(j)*sx(j)+ry(j)*sy(j))*(Srs+Ssr)+(sx(j)^2+sy(j)^2)*Sss);
Kr(:,:,j)=K;

b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M*f(c4n2(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),Kr(:));

A=A+B;

fns = setdiff(1:size(ind4e(:),1), inddb);
u(inddb) = u_D(c4n2(inddb,:));
u(fns) = A(fns,fns)\b(fns);
plot3(c4n2(:,1),c4n2(:,2),u)