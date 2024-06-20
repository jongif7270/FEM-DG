function [error,rate]=FEMHDGerror(t,N)

%%
%f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
%u=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
%%
% ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2));
% uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));
%%
ux=@(x) cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8;
uy=@(x) cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4;
%%
xl=-1;xr=1;yl=-1;yr=1;
% xl=-4;xr=4;yl=-4;yr=4;

iter=t;
% time=1:iter;
error=zeros(1,iter);
Ma=2.^(1:iter);
h=1./Ma;
for i =1:iter
    [~,n4e,ind4e,~,~,~,~,~] = mesh_FEMDG(xl,xr,yl,yr,Ma(i),Ma(i),N);
    [u,V,Dr,Ds,c4n] = FEMHDG2(Ma(i),N);
    error(i)=compute_error_FEMDG(c4n,n4e,ind4e,V,Dr,Ds,u,ux,uy,N);
end

%plot(time,error)
% loglog(2.^time,error)
rate=(log(error(2:iter))-log(error(1:(iter-1))))./(log(h(2:iter))-log(h(1:(iter-1))));