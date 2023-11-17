function [V,M,Dr,Ds,Srr,Srs,Ssr,Sss,c4n,ind4e,inddb,u,A,b,fns,ux,uy] = FEMPOISSON(xl,xr,yl,yr,Mx,My,N)
[c4n,n4e,ind4e,inddb] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V);
number_of_nodes = size(c4n,1);
A = sparse(number_of_nodes,number_of_nodes);
b = zeros(number_of_nodes,1);
u = b;
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;
ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));
I = eye((N+1)*(N+2)/2);
for j=1:size(n4e,1)
xr = (c4n(n4e(j,1),1)-c4n(n4e(j,3),1))/2;
yr = (c4n(n4e(j,1),2)-c4n(n4e(j,3),2))/2;
xs = (c4n(n4e(j,2),1)-c4n(n4e(j,3),1))/2;
ys = (c4n(n4e(j,2),2)-c4n(n4e(j,3),2))/2;
J = xr*ys-xs*yr;
rx=ys/J; ry=-xs/J; sx=-yr/J; sy=xr/J;

M=J*I/(V*V');
Srr=J*(V\Dr)'*(V\Dr);
Srs=J*(V\Dr)'*(V\Ds);
Ssr=J*(V\Ds)'*(V\Dr);
Sss=J*(V\Ds)'*(V\Ds);

A(ind4e(j,:),ind4e(j,:)) = A(ind4e(j,:),ind4e(j,:)) ...
 + J*((rx^2+ry^2)*Srr+(rx*sx+ry*sy)*(Srs+Ssr)+(sx^2+sy^2)*Sss);
b(ind4e(j,:)) = b(ind4e(j,:)) + J*M*f(c4n(ind4e(j,:),:));
end
fns = setdiff(1:number_of_nodes, inddb);
u(inddb) = u_D(c4n(inddb,:));
u(fns) = A(fns,fns)\b(fns);
plot3(c4n(:,1),c4n(:,2),u)
end