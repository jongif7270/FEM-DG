function [u,V2D,Dr,Ds,c4n] = FEMHDG4(M,N)

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

e=1;S=1;k=4*N^2;  f=@(x) e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

e=1;S=1;k=4*N^2;  f=@(x) e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

[~,n4e,ind4e,~,ind4s,e4s,s4e,~,en,ind4p,s4p,e4p,c4n] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

[r1D] = JacobiGL(0,0,N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

[x,y]=Nodes2D(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

xr = (c4n(n4e(:,2),1)-c4n(n4e(:,1),1))/2;
yr = (c4n(n4e(:,2),2)-c4n(n4e(:,1),2))/2;
xs = (c4n(n4e(:,3),1)-c4n(n4e(:,1),1))/2;
ys = (c4n(n4e(:,3),2)-c4n(n4e(:,1),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

h = vecnorm(c4n(n4e(:,3),:) - c4n(n4e(:,2),:),2,2);

n = ((c4n(n4e(:,3),:) - c4n(n4e(:,2),:))./h)*[0 -1 ; 1 0];

[GUGV,GUV,UGV,UV] = compute_matrix(rx,ry,sx,sy,J,h,n,Dr,Ds,M2D,M);




