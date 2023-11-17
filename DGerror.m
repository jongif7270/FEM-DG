function [time,error]=DGerror(xl,xr,yl,yr,t)
iter=t;
time=1:iter;
error=zeros(1,iter);
Ma=2.^(1:iter);
for i =1:iter
    [c4n,n4e,ind4e,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Ma(i),Ma(i),1);
    [~,V,Dr,Ds,u,ux,uy] = DG(xl,xr,yl,yr,Ma(i),Ma(i),1);
    error(i)=computeDGerror(c4n,n4e,ind4e,Dr,Ds,u,ux,uy,V,1);
end
plot(time,error)