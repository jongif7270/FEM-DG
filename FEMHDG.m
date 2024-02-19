function [u,V2D,Dr,Ds,c4n2] = FEMHDG(M,N)

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    
a=@(x) [0,0];b=0;e=1;S=1;k=4*N^2;  
f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +0*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) -0*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); 
ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

