function [A,V2D,Dr,Ds,u,c4n2] = DG2(xl,xr,yl,yr,Mx,My,N)
[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,inddb,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

b = zeros(size(ind4e(:),1),1);
u = b;
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));
I2D = eye((N+1)*(N+2)/2);

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);

[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);

%A = sparse(Mx*My*2*(N+1)*(N+2)/2,Mx*My*2*(N+1)*(N+2)/2);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

Kr=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
    M=J(j)*I2D/(V2D*V2D');
    Srr=J(j)*(V2D\Dr)'*(V2D\Dr);
    Srs=J(j)*(V2D\Dr)'*(V2D\Ds);
    Ssr=J(j)*(V2D\Ds)'*(V2D\Dr);
    Sss=J(j)*(V2D\Ds)'*(V2D\Ds);

    K=J(j)*((rx(j)^2+ry(j)^2)*Srr+(rx(j)*sx(j)+ry(j)*sy(j))*(Srs+Ssr)+(sx(j)^2+sy(j)^2)*Sss);
    Kr(:,:,j)=K;

    b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M*f(c4n2(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),Kr(:));

en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;
sigma=4*(N^2);
for j=1:size(e4s,1)
    if e4s(j,2)~=0
        n=normal(c4n2(ind4s(j,1,1),:)-c4n2(ind4s(j,N+1,1),:));
        h=norm(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
        M=(h/2)*I1D/(V1D*V1D');

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
        
    end
end

fns = setdiff(1:size(ind4e(:),1), inddb);
%fns = 1:size(ind4e(:),1);
u(inddb) = u_D(c4n2(inddb,:));
u(fns) = A(fns,fns)\b(fns);
plot3(c4n2(:,1),c4n2(:,2),u,'.')