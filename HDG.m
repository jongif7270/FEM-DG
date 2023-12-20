function [A,B,C,D,E,b,u] = HDG(M,N)

xl=0;xr=1;yl=0;yr=1;Mx=M;My=M;

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,~,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

b = zeros(size(ind4e(:),1),1);
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);

[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

Kr=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));

for j=1:size(n4e,1)
    M=I2D/(V2D*V2D');
    Srr=(V2D\Dr)'*(V2D\Dr);
    Srs=(V2D\Dr)'*(V2D\Ds);
    Ssr=(V2D\Ds)'*(V2D\Dr);
    Sss=(V2D\Ds)'*(V2D\Ds);

    K=J(j)*((rx(j)^2+ry(j)^2)*Srr+(rx(j)*sx(j)+ry(j)*sy(j))*(Srs+Ssr)+(sx(j)^2+sy(j)^2)*Sss);
    Kr(:,:,j)=K;

    b(ind4e(j,:)) = b(ind4e(j,:)) + J(j)*M*f(c4n2(ind4e(j,:),:));
end

ind=ind4e';
Ir=repmat(ind4e,1,(N+1)*(N+2)/2)';
Jr=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ir(:),Jr(:),Kr(:));

B=zeros(Mx*My*(N+1)*(N+2),size(e4s,1)*(N+1));

C=B';

D=zeros(size(e4s,1)*(N+1),size(e4s,1)*(N+1));

en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;
k=4*N^2;

for j=1:size(e4s,1)
    if e4s(j,2)~=0
        n=normal(c4n2(ind4s(j,N+1,1),:)-c4n2(ind4s(j,1,1),:));
        h=norm(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
        M=I1D/(V1D*V1D');

        A(ind4e(e4s(j,1),:),ind4s(j,:,1))=A(ind4e(e4s(j,1),:),ind4s(j,:,1))-(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        A(ind4e(e4s(j,2),:),ind4s(j,:,2))=A(ind4e(e4s(j,2),:),ind4s(j,:,2))+(h/2)*((rx(e4s(j,2))*Dr(en(j,:,2),:)'*M+sx(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(1)+(ry(e4s(j,2))*Dr(en(j,:,2),:)'*M+sy(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(2));
        B(ind4e(e4s(j,1),:),(j-1)*(N+1)+(1:N+1))=B(ind4e(e4s(j,1),:),(j-1)*(N+1)+(1:N+1))+(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        B(ind4e(e4s(j,2),:),(j-1)*(N+1)+(1:N+1))=B(ind4e(e4s(j,2),:),(j-1)*(N+1)+(1:N+1))-(h/2)*((rx(e4s(j,2))*Dr(en(j,:,2),:)'*M+sx(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(1)+(ry(e4s(j,2))*Dr(en(j,:,2),:)'*M+sy(e4s(j,2))*Ds(en(j,:,2),:)'*M)*n(2));
        
        A(ind4s(j,:,1),ind4e(e4s(j,1),:))=A(ind4s(j,:,1),ind4e(e4s(j,1),:))-(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        A(ind4s(j,:,2),ind4e(e4s(j,2),:))=A(ind4s(j,:,1),ind4e(e4s(j,2),:))+(h/2)*((rx(e4s(j,2))*M*Dr(en(j,:,2),:)+sx(e4s(j,2))*M*Ds(en(j,:,2),:))*n(1)+(ry(e4s(j,2))*M*Dr(en(j,:,2),:)+sy(e4s(j,2))*M*Ds(en(j,:,2),:))*n(2));
        C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,1),:))=C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,1),:))+(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,2),:))=C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,2),:))-(h/2)*((rx(e4s(j,2))*M*Dr(en(j,:,2),:)+sx(e4s(j,2))*M*Ds(en(j,:,2),:))*n(1)+(ry(e4s(j,2))*M*Dr(en(j,:,2),:)+sy(e4s(j,2))*M*Ds(en(j,:,2),:))*n(2));
    else
        n=normal(c4n2(ind4s(j,N+1,1),:)-c4n2(ind4s(j,1,1),:));
        h=norm(c4n(n4s(j,2),:)-c4n(n4s(j,1),:));
        M=I1D/(V1D*V1D');

        A(ind4e(e4s(j,1),:),ind4s(j,:,1))=A(ind4e(e4s(j,1),:),ind4s(j,:,1))-(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        B(ind4e(e4s(j,1),:),(j-1)*(N+1)+(1:N+1))=B(ind4e(e4s(j,1),:),(j-1)*(N+1)+(1:N+1))+(h/2)*((rx(e4s(j,1))*Dr(en(j,:,1),:)'*M+sx(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(1)+(ry(e4s(j,1))*Dr(en(j,:,1),:)'*M+sy(e4s(j,1))*Ds(en(j,:,1),:)'*M)*n(2));
        
        A(ind4s(j,:,1),ind4e(e4s(j,1),:))=A(ind4s(j,:,1),ind4e(e4s(j,1),:))-(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
        C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,1),:))=C((j-1)*(N+1)+(1:N+1),ind4e(e4s(j,1),:))+(h/2)*((rx(e4s(j,1))*M*Dr(en(j,:,1),:)+sx(e4s(j,1))*M*Ds(en(j,:,1),:))*n(1)+(ry(e4s(j,1))*M*Dr(en(j,:,1),:)+sy(e4s(j,1))*M*Ds(en(j,:,1),:))*n(2));
    end
end

for j=1:size(n4e,1)
    he=norm(c4n2(ind4s(s4e(j,1),N+1,1),:)-c4n2(ind4s(s4e(j,1),1,1),:));
    for i=1:3
        h=norm(c4n2(ind4s(s4e(j,i),N+1,1),:)-c4n2(ind4s(s4e(j,i),1,1),:));
        M=I1D/(V1D*V1D');

        A(ind4s(s4e(j,i),:,1),ind4s(s4e(j,i),:,1))=A(ind4s(s4e(j,i),:,1),ind4s(s4e(j,i),:,1))+(h/2)*(k/he)*M;

        B(ind4s(s4e(j,i),:,1),(s4e(j,i)-1)*(N+1)+(1:N+1))=B(ind4s(s4e(j,i),:,1),(s4e(j,i)-1)*(N+1)+(1:N+1))-(h/2)*(k/he)*M;

        C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4s(s4e(j,i),:,1))=C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4s(s4e(j,i),:,1))-(h/2)*(k/he)*M;

        D((s4e(j,i)-1)*(N+1)+(1:N+1),(s4e(j,i)-1)*(N+1)+(1:N+1))=D((s4e(j,i)-1)*(N+1)+(1:N+1),(s4e(j,i)-1)*(N+1)+(1:N+1))+(h/2)*(k/he)*M;
    end
end


E=[A B; C D];

V=-(D-(C/A)*B)\(C/A)*b;
u=A\(b-B*V);

spy(E);

%spy(A)
% fns = setdiff(1:size(ind4e(:),1), inddb);
% u(inddb) = u_D(c4n2(inddb,:));
%fns = 1:size(ind4e(:),1);
% u(fns) = A(fns,fns)\b(fns);
%u=A\b;
%plot3(c4n2(:,1),c4n2(:,2),u,'.')