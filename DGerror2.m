function [time,error]=DGerror2(xl,xr,yl,yr,M,t)

%%%%
%f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
%u=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));
%%%%

iter=t;
time=1:iter;
error=zeros(1,iter);
Na=2.^(1:iter);
for i =1:iter
    [c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,M,M,Na(i));
    [ind4e,~,c4n2] = indexforDG2(xl,xr,yl,yr,M,M,Na(i));
    [~,V,Dr,Ds,u,~] = DG3(xl,xr,yl,yr,M,M,Na(i));
    error(i)=computeDGerror(c4n,c4n2,n4e,ind4e,Dr,Ds,u,ux,uy,V,Na(i));
end
plot(time,error)