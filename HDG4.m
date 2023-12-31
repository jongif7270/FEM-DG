function [u,V2D,Dr,Ds,c4n2] = HDG4(M,N)

%% 
%xl=0;xr=1;yl=0;yr=1;Mx=M;My=M;    a=[0,0];b=0;e=1;S=1;k=4*N^2;  f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));

%% Ex 3.3
xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=[0.8,0.6];b=1;e=0;S=1;k=4*N^2;  f=@(x) sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4)); %u=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)

%% Ex 3.4
%xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [exp(x(:,1))*(x(:,2)*cos(x(:,2))+sin(x(:,2))),-exp(x(:,1))*x(:,2)*sin(x(:,2))];b=0;e=1;S=-1;k=1;

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,~,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

d = zeros(size(ind4e(:),1),1);

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
s4e = computeS4e(n4e);

count=sum(e4s(:, 2) == 0);
inddb2=zeros(1,count*(N+1));
tmp=0;
for i=1:size(e4s,1)
    if e4s(i,2)==0
        tmp=tmp+1;
        inddb2((tmp-1)*(N+1)+(1:N+1))=(i-1)*(N+1)+(1:N+1);
    end
end

inddb3=[];
inddb4=[];
for i=1:size(e4s,1)
    if e4s(i,2)==0
        if mod(e4s(i,1),2)~=0
            inddb3=[inddb3 (i-1)*(N+1)+(1:N+1)];
        else
            inddb4=[inddb4 (i-1)*(N+1)+(1:N+1)];
        end
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

Aa=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2,size(n4e,1));
Ab=Aa;
Ac=Aa;
Ad=Aa;
Ar=Aa;

Ba=zeros((N+1)*(N+2)/2,3*(N+1),size(n4e,1));
Bb=zeros((N+1)*(N+2)/2,3*(N+1),size(n4e,1));
Br=zeros((N+1)*(N+2)/2,3*(N+1),size(n4e,1));

Ca=zeros(3*(N+1),(N+1)*(N+2)/2,size(n4e,1));
Cb=zeros(3*(N+1),(N+1)*(N+2)/2,size(n4e,1));
Cr=zeros(3*(N+1),(N+1)*(N+2)/2,size(n4e,1));

Da=zeros(3*(N+1),3*(N+1),size(n4e,1));

Fr=zeros(3*(N+1),3*(N+1),size(n4e,1));

Gr=zeros(3*(N+1),1,size(n4e,1));

T=[(s4e(:,1)-1)*(N+1)+(1:N+1) (s4e(:,2)-1)*(N+1)+(1:N+1) (s4e(:,3)-1)*(N+1)+(1:N+1)];

en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

fns = setdiff(1:size(e4s,1)*(N+1), inddb2);
%fns = 1:size(e4s,1)*(N+1);
V(inddb3) = zeros(length(inddb3),1);
V(inddb4) = zeros(length(inddb4),1);

for j=1:size(n4e,1)
    Aa(:,:,j)=e*J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds)+...
              b*J(j)*M2D+...
              -J(j)*(a(1)*(rx(j)*Dr'+sx(j)*Ds')*M2D+a(2)*(ry(j)*Dr'+sy(j)*Ds')*M2D);
    da=J(j)*M2D*f(c4n2(ind4e(j,:),:));
    ht=norm(c4n(n4e(j,2),:)-c4n(n4e(j,1),:));
    for i=1:size(n4e,2)
        n=normal(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));
        h=norm(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));

        if dot(a,n)>0
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1),j)=e*(h/2)*(k/ht)*M;
    
            Ab(:,en(s4e(j,i),:,1),j)=Ab(:,en(s4e(j,i),:,1),j)-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:,j)=Ac(en(s4e(j,i),:,1),:,j)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1),j)=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1),j)+e*(h/2)*(k/ht)*M+dot(a,n)*(h/2)*M;
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1),j)=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1),j)=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1),j)-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:,j)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1),j)=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1),j)-e*(h/2)*(k/ht)*M-dot(a,n)*(h/2)*M;
            else
                Ba(:,(i-1)*(N+1)+(1:N+1),j)=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1),j)=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1),j)-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(1:N+1),:,j)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1),j)=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1),j)-e*(h/2)*(k/ht)*M-dot(a,n)*(h/2)*M;
            end
        else
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1),j)=e*(h/2)*(k/ht)*M-dot(a,n)*(h/2)*M;
    
            Ab(:,en(s4e(j,i),:,1),j)=Ab(:,en(s4e(j,i),:,1),j)-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:,j)=Ac(en(s4e(j,i),:,1),:,j)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1),j)=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1),j)+e*(h/2)*(k/ht)*M;
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1),j)=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1),j)=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1),j)-e*(h/2)*(k/ht)*M+dot(a,n)*(h/2)*M;
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:,j)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1),j)=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1),j)-e*(h/2)*(k/ht)*M;
            else
                Ba(:,(i-1)*(N+1)+(1:N+1),j)=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1),j)=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1),j)-e*(h/2)*(k/ht)*M+dot(a,n)*(h/2)*M;
    
                Ca((i-1)*(N+1)+(1:N+1),:,j)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1),j)=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1),j)-e*(h/2)*(k/ht)*M;
            end
        end
    end
    Ar(:,:,j)=Aa(:,:,j)+Ab(:,:,j)+Ac(:,:,j)+Ad(:,:,j);
    Br(:,:,j)=Ba(:,:,j)+Bb(:,:,j);
    Cr(:,:,j)=Ca(:,:,j)+Cb(:,:,j);

    Fr(:,:,j)=Cr(:,:,j)*(Ar(:,:,j)\Br(:,:,j))-Da(:,:,j);
    Gr(:,1,j)=Cr(:,:,j)*(Ar(:,:,j)\da);

    d(ind4e(j,:))=d(ind4e(j,:))+da;
end

ind=ind4e';
TA=T';

Ia=repmat(ind4e,1,(N+1)*(N+2)/2)';
Ja=(repmat(ind(:),1,(N+1)*(N+2)/2))';
A=sparse(Ia(:),Ja(:),Ar(:));

Ib=repmat(ind4e,1,3*(N+1))';
Jb=(repmat(TA(:),1,(N+1)*(N+2)/2))';
B=sparse(Ib(:),Jb(:),Br(:));

%Ic=(repmat(T,1,(N+1)*(N+2)/2))';
%Jc=repmat(ind(:),1,3*(N+1))';
%C=sparse(Ic(:),Jc(:),Cr(:));

%Id=repmat(T,1,3*(N+1))';
%Jd=(repmat(TA(:),1,3*(N+1)))';
%D=sparse(Id(:),Jd(:),Da(:));

If=repmat(T,1,3*(N+1))';
Jf=(repmat(TA(:),1,3*(N+1)))';
F=sparse(If(:),Jf(:),Fr(:));

Ig=TA;
Jg=ones(3*(N+1),size(s4e,1));
G=sparse(Ig(:),Jg(:),Gr(:));

%F=C*(A\B)-D;
%G=C*(A\d);
V(fns)=F(fns,fns)\G(fns);

u=A\(d-(B*V'));

plot3(c4n2(:,1),c4n2(:,2),u,'.')