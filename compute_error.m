function [time,error]=compute_error(k,t)

xl=0;xr=1;yl=0;yr=1;

ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));

iter=t;
time=1:iter;
error=zeros(1,iter);
Ma=2.^(1:iter);
for i =1:iter
    [c4n,n4e,ind4e,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Ma(i),Ma(i),k);
    [~,u,Dr,Ds,V] = FEM(xl,xr,yl,yr,Ma(i),Ma(i),k);
    error(i)=compute_error_fem_2d_triangle(c4n,n4e,ind4e,V,Dr,Ds,u,ux,uy,k);
end
plot(time,error)