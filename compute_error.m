function [time,error]=compute_error(xl,xr,yl,yr,k,t)
iter=t;
time=1:iter;
error=zeros(1,iter);
Ma=2.^(1:iter);
for i =1:iter
    [c4n,n4e,ind4e,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Ma(i),Ma(i),k);
    [~,M,Dr,Ds,~,~,~,~,~,~,~,u,~,~,~,ux,uy] = FEMPOISSON(xl,xr,yl,yr,Ma(i),Ma(i),k);
    error(i)=compute_error_fem_2d_triangle2(c4n,n4e,ind4e,M,Dr,Ds,u,ux,uy);
end
plot(time,error)