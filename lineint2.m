function [A]=lineint2(xl,xr,yl,yr,Mx,My,N)
[r] = Nodes1D_equi(N);
[V] = Vandermonde1D(N,r);
I = eye(size(V,1));

[ind4e,~,~,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);

[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);

%A = zeros(Mx*My*2*(N+1)*(N+2)/2,Mx*My*2*(N+1)*(N+2)/2);
A = sparse(Mx*My*2*(N+1)*(N+2)/2,Mx*My*2*(N+1)*(N+2)/2);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;
sigma=4*(N^2);
for j=1:size(e4s,1)
    if e4s(j,2)~=0
    n=normal(c4n2(ind4s(j,1,1),:)-c4n2(ind4s(j,N+1,1),:));
    %n=normal(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
    h=norm(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
    M=(h/2)*I/(V*V');

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
    end
end