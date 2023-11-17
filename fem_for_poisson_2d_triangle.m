function [u,A,b,fns] = fem_for_poisson_2d_triangle(c4n,n4e,n4db, ...
ind4e,M_R,Srr_R,Srs_R,Ssr_R,Sss_R,f,u_D)
number_of_nodes = size(c4n,1);
A = sparse(number_of_nodes,number_of_nodes);
b = zeros(number_of_nodes,1);
u = b;
for j=1:size(n4e,1)
xr = (c4n(n4e(j,1),1)-c4n(n4e(j,3),1))/2;
yr = (c4n(n4e(j,1),2)-c4n(n4e(j,3),2))/2;
xs = (c4n(n4e(j,2),1)-c4n(n4e(j,3),1))/2;
ys = (c4n(n4e(j,2),2)-c4n(n4e(j,3),2))/2;
J = xr*ys-xs*yr;
rx=ys/J; ry=-xs/J; sx=-yr/J; sy=xr/J;

A(ind4e(j,:),ind4e(j,:)) = A(ind4e(j,:),ind4e(j,:)) ...
+J*((rx^2+ry^2)*Srr_R+(rx*sx+ry*sy)*(Srs_R+Ssr_R)+(sx^2+sy^2)*Sss_R);

b(ind4e(j,:)) = b(ind4e(j,:)) + J*M_R*f(c4n(ind4e(j,:),:));
end
fns = setdiff(1:number_of_nodes, n4db);
u(n4db) = u_D(c4n(n4db,:));
u(fns) = A(fns,fns)\b(fns);
end
