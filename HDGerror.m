function [time,error]=HDGerror(t,N)

%%%%
%f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
%u=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));
%%%%

xl=0;xr=1;yl=0;yr=1;

iter=t;
time=1:iter;
error=zeros(1,iter);
Ma=2.^(1:iter);
for i =1:iter
    [c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Ma(i),Ma(i),N);
    [ind4e,~,c4n2] = indexforDG2(xl,xr,yl,yr,Ma(i),Ma(i),N);
    [u,V,Dr,Ds,~] = HDG4(Ma(i),N);
    error(i)=computeDGerror(c4n,c4n2,n4e,ind4e,Dr,Ds,u,ux,uy,V,N);
end
%plot(time,error)
loglog(2.^time,error)