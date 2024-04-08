function [u,V2D,Dr,Ds,c4n] = FEMDG(M,N)

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); 
ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;

[c4n,n4e,ind4e,inddb,ind4s,e4s] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

fns=setdiff(1:size(c4n),inddb);

xr = (c4n(n4e(:,2),1)-c4n(n4e(:,1),1))/2;
yr = (c4n(n4e(:,2),2)-c4n(n4e(:,1),2))/2;
xs = (c4n(n4e(:,3),1)-c4n(n4e(:,1),1))/2;
ys = (c4n(n4e(:,3),2)-c4n(n4e(:,1),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

b = zeros(size(c4n,1),1);

S=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
    S(:,:,j)=J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds);
    b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M2D*f(c4n(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),S(:));

u(inddb) = u_D(c4n(inddb,:));
u(fns) = A(fns,fns)\b(fns);

%%
T = delaunay(c4n(:,1),c4n(:,2));
trisurf(T,c4n(:,1),c4n(:,2),u,"Linestyle","none")
% plot3(c4n(:,1),c4n(:,2),u,'.')