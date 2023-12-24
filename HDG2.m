function u = HDG2(M,N)

xl=0;xr=1;yl=0;yr=1;Mx=M;My=M;

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,inddb,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

b = zeros(size(ind4e(:),1),1);
f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
u_D=@(x) x(:,1)*0;

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);

inddb2=[];
for i=1:size(e4s,1)
    if e4s(i,2)==0
    inddb2=[inddb2 (i-1)*(N+1)+(1:N+1)];
    end
end

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

M=I1D/(V1D*V1D');

for j=1:size(n4e,1)
    ht=norm(c4n2(ind4s(s4e(j,1),N+1,1),:)-c4n2(ind4s(s4e(j,1),1,1),:));
    for i=1:size(n4e,2)
        n=normal(c4n2(ind4s(s4e(j,i),N+1,1),:)-c4n2(ind4s(s4e(j,i),1,1),:));
        h=norm(c4n2(ind4s(s4e(j,i),N+1,1),:)-c4n2(ind4s(s4e(j,i),1,1),:));

        D((s4e(j,i)-1)*(N+1)+(1:N+1),(s4e(j,i)-1)*(N+1)+(1:N+1))=D((s4e(j,i)-1)*(N+1)+(1:N+1),(s4e(j,i)-1)*(N+1)+(1:N+1))+(h/2)*(k/ht)*M;
        
        if e4s(s4e(j,i),2)==j
            n=-n;
            A(ind4e(j,:),fliplr(ind4s(s4e(j,i),:,2)))=A(ind4e(j,:),fliplr(ind4s(s4e(j,i),:,2)))-(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            A(fliplr(ind4s(s4e(j,i),:,2)),ind4e(j,:))=A(fliplr(ind4s(s4e(j,i),:,2)),ind4e(j,:))-(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            A(fliplr(ind4s(s4e(j,i),:,2)),fliplr(ind4s(s4e(j,i),:,2)))=A(fliplr(ind4s(s4e(j,i),:,2)),fliplr(ind4s(s4e(j,i),:,2)))+(h/2)*(k/ht)*M;

            B(ind4e(j,:),(s4e(j,i)-1)*(N+1)+(N+1:-1:1))=B(ind4e(j,:),(s4e(j,i)-1)*(N+1)+(N+1:-1:1))+(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            B(fliplr(ind4s(s4e(j,i),:,2)),(s4e(j,i)-1)*(N+1)+(N+1:-1:1))=B(fliplr(ind4s(s4e(j,i),:,2)),(s4e(j,i)-1)*(N+1)+(N+1:-1:1))-(h/2)*(k/ht)*M;

            C((s4e(j,i)-1)*(N+1)+(N+1:-1:1),ind4e(j,:))=C((s4e(j,i)-1)*(N+1)+(N+1:-1:1),ind4e(j,:))+(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            C((s4e(j,i)-1)*(N+1)+(N+1:-1:1),fliplr(ind4s(s4e(j,i),:,2)))=C((s4e(j,i)-1)*(N+1)+(N+1:-1:1),fliplr(ind4s(s4e(j,i),:,2)))-(h/2)*(k/ht)*M;
        else
            A(ind4e(j,:),ind4s(s4e(j,i),:,1))=A(ind4e(j,:),ind4s(s4e(j,i),:,1))-(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            A(ind4s(s4e(j,i),:,1),ind4e(j,:))=A(ind4s(s4e(j,i),:,1),ind4e(j,:))-(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            A(ind4s(s4e(j,i),:,1),ind4s(s4e(j,i),:,1))=A(ind4s(s4e(j,i),:,1),ind4s(s4e(j,i),:,1))+(h/2)*(k/ht)*M;
        
            B(ind4e(j,:),(s4e(j,i)-1)*(N+1)+(1:N+1))=B(ind4e(j,:),(s4e(j,i)-1)*(N+1)+(1:N+1))+(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            B(ind4s(s4e(j,i),:,1),(s4e(j,i)-1)*(N+1)+(1:N+1))=B(ind4s(s4e(j,i),:,1),(s4e(j,i)-1)*(N+1)+(1:N+1))-(h/2)*(k/ht)*M;

            C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4e(j,:))=C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4e(j,:))+(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4s(s4e(j,i),:,1))=C((s4e(j,i)-1)*(N+1)+(1:N+1),ind4s(s4e(j,i),:,1))-(h/2)*(k/ht)*M;
        end
    end
end


E=[A B; C D];

spy(E);
fns = setdiff(1:size(e4s,1)*(N+1), inddb2);
V(inddb2) = zeros(length(inddb2),1);
F=C*(A\B)-D;
G=C*(A\b);
V(fns)=F(fns,fns)\G(fns);

u=A\(b-(B*V'));

%u(inddb) = u_D(c4n2(inddb,:));

%Iu=eye(Mx*My*(N+1)*(N+2));
%u=(Iu-A\B/D*C)\A\b;

plot3(c4n2(:,1),c4n2(:,2),u,'.')