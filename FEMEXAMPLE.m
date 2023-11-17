function[error,time]=FEMEXAMPLE
iter = 6;
k = 1;
xl = 0; xr = 1; yl = 0; yr = 1; M = 2.^(1:iter);
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;
ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));

error=zeros(1,iter);
time=zeros(1,iter);
h=1./M;
for j=1:iter
[c4n,n4e,ind4e,n4db] = mesh_fem_2d_triangle(xl,xr,yl,yr,M(j),M(j),k);
[M_R,Srr_R,Srs_R,Ssr_R,Sss_R,Dr_R,Ds_R] = get_matrices_2d_triangle(k);
u = fem_for_poisson_2d_triangle(c4n,n4e,n4db,ind4e,M_R,Srr_R,Srs_R, ...
Ssr_R,Sss_R,f,u_D);
error(j) = compute_error_fem_2d_triangle2(c4n,n4e,ind4e,M_R,Dr_R, ...
Ds_R,u,ux,uy);
end
