function [u,V2D,Dr,Ds,c4n] = FEMDG(M,N)

Mx=M;My=M; 

f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% f=@(x) (sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);


xl=-1;xr=1;yl=-1;yr=1;

% xl=-4;xr=4;yl=-4;yr=4;

[~,n4e,ind4e,inddb,ind4s,e4s,~,n4s,en,~,~,~,c4n] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

[r1D] = JacobiGL(0,0,N);
% [r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

[x,y]=Nodes2D(N);
% [x,y]=Nodes2D_equi(N);
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
u=b;

S=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
    S(:,:,j)=J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds);
    b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M2D*f(c4n(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),S(:));

sigma=4*N^2;
for j=1:size(e4s,1)
    n=normal(c4n(ind4s(j,end,1),:)-c4n(ind4s(j,1,1),:));
    h=norm(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
    if e4s(j,2)~=0
        A(ind4s(j,:,1),ind4e(e4s(j,1),:))=A(ind4s(j,:,1),ind4e(e4s(j,1),:))+(-1/2)*(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        A(ind4s(j,:,2),ind4e(e4s(j,1),:))=A(ind4s(j,:,2),ind4e(e4s(j,1),:))+(1/2)*(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        A(ind4s(j,:,1),ind4e(e4s(j,2),:))=A(ind4s(j,:,1),ind4e(e4s(j,2),:))+(-1/2)*(h/2)*((rx(e4s(j,2))*M*Dr(en(j,:,2),:)+sx(e4s(j,2))*M*Ds(en(j,:,2),:))*n(1)+(ry(e4s(j,2))*M*Dr(en(j,:,2),:)+sy(e4s(j,2))*M*Ds(en(j,:,2),:))*n(2));
        A(ind4s(j,:,2),ind4e(e4s(j,2),:))=A(ind4s(j,:,2),ind4e(e4s(j,2),:))+(1/2)*(h/2)*((rx(e4s(j,2))*M*Dr(en(j,:,2),:)+sx(e4s(j,2))*M*Ds(en(j,:,2),:))*n(1)+(ry(e4s(j,2))*M*Dr(en(j,:,2),:)+sy(e4s(j,2))*M*Ds(en(j,:,2),:))*n(2));

        A(ind4e(e4s(j,1),:),ind4s(j,:,1))=A(ind4e(e4s(j,1),:),ind4s(j,:,1))+(-1/2)*(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        A(ind4e(e4s(j,2),:),ind4s(j,:,1))=A(ind4e(e4s(j,2),:),ind4s(j,:,1))+(-1/2)*(h/2)*((rx(e4s(j,2))*Dr(en(j,:,2),:)'*M+sx(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(1)+(ry(e4s(j,2))*Dr(en(j,:,2),:)'*M+sy(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(2));
        A(ind4e(e4s(j,1),:),ind4s(j,:,2))=A(ind4e(e4s(j,1),:),ind4s(j,:,2))+(1/2)*(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        A(ind4e(e4s(j,2),:),ind4s(j,:,2))=A(ind4e(e4s(j,2),:),ind4s(j,:,2))+(1/2)*(h/2)*((rx(e4s(j,2))*Dr(en(j,:,2),:)'*M+sx(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(1)+(ry(e4s(j,2))*Dr(en(j,:,2),:)'*M+sy(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(2));

        A(ind4s(j,:,1),ind4s(j,:,1))=A(ind4s(j,:,1),ind4s(j,:,1))+(sigma/h)*(h/2)*M;
        A(ind4s(j,:,2),ind4s(j,:,1))=A(ind4s(j,:,2),ind4s(j,:,1))-(sigma/h)*(h/2)*M;
        A(ind4s(j,:,1),ind4s(j,:,2))=A(ind4s(j,:,1),ind4s(j,:,2))-(sigma/h)*(h/2)*M;
        A(ind4s(j,:,2),ind4s(j,:,2))=A(ind4s(j,:,2),ind4s(j,:,2))+(sigma/h)*(h/2)*M;
    else
        A(ind4s(j,:,1),ind4e(e4s(j,1),:))=A(ind4s(j,:,1),ind4e(e4s(j,1),:))-(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        A(ind4e(e4s(j,1),:),ind4s(j,:,1))=A(ind4e(e4s(j,1),:),ind4s(j,:,1))-(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));

        b(ind4e(e4s(j,1),:))=b(ind4e(e4s(j,1),:))-(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2))*ue(c4n(ind4s(j,:,1),:));
    end
end

% u(inddb) = ue(c4n(inddb,:));
% u(fns) = A(fns,fns)\b(fns);
u = A\b;

% spy(A)

%%
T = delaunay(c4n(:,1),c4n(:,2));
trisurf(T,c4n(:,1),c4n(:,2),u,"Linestyle","none")
% plot3(c4n(:,1),c4n(:,2),u,'.')