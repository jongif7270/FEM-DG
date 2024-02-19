function [u,V2D,Dr,Ds,c4n2] = HDG7

M=8;N=8;

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))), -exp(x(:,1)).*x(:,2).*sin(x(:,2))]; b=0;e=1;S=1;k=4*N^2; f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +(exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2)))).*pi.*cos(pi*x(:,1)).*sin(pi*x(:,2)) +(-exp(x(:,1)).*x(:,2).*sin(x(:,2))).*pi.*sin(pi*x(:,1)).*cos(pi*x(:,2)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,~,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);
[~,~,~,c4na,~] = indexforDG(xl,xr,yl,yr,Mx,My,2*N);
ind4sa=edge(n4e,2*N);

d = zeros(size(ind4e(:),1),1);

[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
[V1DA] = Vandermonde1D(2*N,r1D);
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

inddb2=zeros(size(e4s,1)*(N+1));

for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        if e4s(s4e(j,i),2)==0
            inddb2((s4e(j,i)-1)*(N+1)+(1:N+1))=(s4e(j,i)-1)*(N+1)+(1:N+1);
        end
    end
end
inddb2=nonzeros(inddb2)';


[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
V2DA=Vandermonde2D(2*N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

Da=zeros(3*(N+1),3*(N+1));

Fr=zeros(3*(N+1),3*(N+1),size(n4e,1));

Gr=zeros(3*(N+1),size(n4e,1));

AFr=zeros((N+1)*(N+2)/2,size(n4e,1));

ABr=zeros((N+1)*(N+2)/2,3*(N+1),size(n4e,1));

V=zeros(size(e4s,1)*(N+1),1);

T=[(s4e(:,1)-1)*(N+1)+(1:N+1) (s4e(:,2)-1)*(N+1)+(1:N+1) (s4e(:,3)-1)*(N+1)+(1:N+1)];

en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

fns = setdiff(1:size(e4s,1)*(N+1), inddb2);

ae=zeros((N+1)*(N+2)/2,2,size(n4e,1));
for j=1:size(n4e,1)
    for l=1:(N+1)*(N+2)/2
        ae(l,:,j)=a(c4n2(((N+1)*(N+2)/2)*(j-1)+l,:));
    end
end

al=zeros(N+1,2,size(e4s,1));
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        for l=1:N+1
            al(l,:,s4e(j,i))=a(c4n2(ind4s(s4e(j,i),l,1),:));
        end
    end
end

aea=zeros((2*N+1)*(2*N+2)/2,2,size(n4e,1));
for j=1:size(n4e,1)
    for l=1:(2*N+1)*(2*N+2)/2
        aea(l,:,j)=a(c4na(((2*N+1)*(2*N+2)/2)*(j-1)+l,:));
    end
end

ala=zeros(2*N+1,2,size(e4s,1));
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        for l=1:2*N+1
            ala(l,:,s4e(j,i))=a(c4na(ind4sa(s4e(j,i),l,1),:));
        end
    end
end

for j=1:size(n4e,1)
    Ab=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);
    Ac=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);
    Ad=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);

    Ba=zeros((N+1)*(N+2)/2,3*(N+1));
    Bb=zeros((N+1)*(N+2)/2,3*(N+1));

    Ca=zeros(3*(N+1),(N+1)*(N+2)/2);
    Cb=zeros(3*(N+1),(N+1)*(N+2)/2);

    Aa=e*J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds)+...
              b*J(j)*M2D+...
              -J(j)*((rx(j)*Dr'+sx(j)*Ds')*M2D.*(V2DA*aea(:,1,j))'+(ry(j)*Dr'+sy(j)*Ds')*M2D.*(V2DA*aea(:,2,j))');
    da=J(j)*M2D*f(c4n2(ind4e(j,:),:));
    ht=norm(c4n(n4e(j,2),:)-c4n(n4e(j,1),:));
    for i=1:size(n4e,2)
        n=normal(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));
        h=norm(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));

        v=ala(:,:,s4e(j,i))*n;
        vp=0;
        vn=0;
        for l=1:2*N+1
            if v(l)>0
                vp=vp+1;
            elseif v(l)<0
                vn=vn+1;
            end
        end

        if vp>vn
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=e*(h/2)*(k/ht)*M;
    
            Ab(:,en(s4e(j,i),:,1))=Ab(:,en(s4e(j,i),:,1))-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:)=Ac(en(s4e(j,i),:,1),:)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))+e*(h/2)*(k/ht)*M+(h/2)*M.*(V1DA*v)';
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M-(h/2)*M.*(V1DA*v)';
            else
                Ba(:,(i-1)*(N+1)+(1:N+1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(1:N+1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M-(h/2)*M.*(V1DA*v)';
            end
            if e4s(s4e(j,i),2)==0
                da=da-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2))*ue(c4n2(ind4s(s4e(j,i),:,1),:));
                da(en(s4e(j,i),:,1))=da(en(s4e(j,i),:,1))-(-e*(h/2)*(k/ht)*M)*ue(c4n2(ind4s(s4e(j,i),:,1),:));
            end
        else
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=e*(h/2)*(k/ht)*M-(h/2)*M.*(V1DA*v)';
    
            Ab(:,en(s4e(j,i),:,1))=Ab(:,en(s4e(j,i),:,1))-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:)=Ac(en(s4e(j,i),:,1),:)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))+e*(h/2)*(k/ht)*M;
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))-e*(h/2)*(k/ht)*M+(h/2)*M.*(V1DA*v)';
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M;
            else
                Ba(:,(i-1)*(N+1)+(1:N+1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))-e*(h/2)*(k/ht)*M+(h/2)*M.*(V1DA*v)';
    
                Ca((i-1)*(N+1)+(1:N+1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M;
            end

            if e4s(s4e(j,i),2)==0
                da=da-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2))*ue(c4n2(ind4s(s4e(j,i),:,1),:));
                da(en(s4e(j,i),:,1))=da(en(s4e(j,i),:,1))-(-e*(h/2)*(k/ht)*M+(h/2)*M.*(V1DA*v)')*ue(c4n2(ind4s(s4e(j,i),:,1),:));
            end
        end
    end
    Ar=Aa+Ab+Ac+Ad;
    Br=Ba+Bb;
    Cr=Ca+Cb;

    Fr(:,:,j)=-Cr*(Ar\Br)+Da;
    Gr(:,j)=-Cr*(Ar\da);

    AFr(:,j)=Ar\da;
    ABr(:,:,j)=Ar\Br;

    d(ind4e(j,:))=d(ind4e(j,:))+da;
end

ind=ind4e';
TA=T';

If=repmat(T,1,3*(N+1))';
Jf=(repmat(TA(:),1,3*(N+1)))';
F=sparse(If(:),Jf(:),Fr(:));

Ig=TA;
Jg=ones(3*(N+1),size(s4e,1));
G=sparse(Ig(:),Jg(:),Gr(:));

V(fns)=F(fns,fns)\G(fns);

Iaf=ind;
Jaf=ones((N+1)*(N+2)/2,size(s4e,1));
AF=sparse(Iaf(:),Jaf(:),AFr(:));

Iab=repmat(ind4e,1,3*(N+1))';
Jab=(repmat(TA(:),1,(N+1)*(N+2)/2))';
AB=sparse(Iab(:),Jab(:),ABr(:));

u=AF-AB*V;

plot3(c4n2(:,1),c4n2(:,2),u,'.')