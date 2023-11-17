function [time,error]=compute_error2(xl,xr,yl,yr,Mx,My,k)
iter=k;
time=1:iter;
error=zeros(1,iter);
Ka=2.^(1:iter);
for i =1:iter
    [c4n,n4e,ind4e,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,Ka(i));
    [~,M,Dr,Ds,~,~,~,~,~,~,~,u,~,~,~,ux,uy] = FEMPOISSON(xl,xr,yl,yr,Mx,My,Ka(i));
    error(i)=compute_error_fem_2d_triangle2(c4n,n4e,ind4e,M,Dr,Ds,u,ux,uy);
end
plot(time,error)