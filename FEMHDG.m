function [u,V2D,Dr,Ds,c4n] = FEMHDG(M,N)


% xl=-4;xr=4;yl=-4;yr=4;Mx=M;My=M;

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

% a=@(x) [0,0];b=0;e=1;S=-1;k=1;  f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +0*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) +0*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

a=@(x) [0.8,0.6];b=1;e=1;S=1;k=4*N^2;  f=@(x) b.*(1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)) +0.8*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8 +0.6*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4 +e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

% a=@(x) [exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))), -exp(x(:,1)).*x(:,2).*sin(x(:,2))]; b=0;e=1;S=-1;k=1; f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +pi*exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))).*cos(pi*x(:,1)).*sin(pi*x(:,2)) +pi*(-exp(x(:,1)).*x(:,2).*sin(x(:,2))).*sin(pi*x(:,1)).*cos(pi*x(:,2)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

[c4n,n4e,ind4e,~,ind4s,e4s,s4e,~,en,~,~,~] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

% [r1D] = JacobiGL(0,0,N);
[r1D] = Nodes1D_equi(N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

% [x,y]=Nodes2D(N);
[x,y]=Nodes2D_equi(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

inddb2=zeros(size(e4s,1)*(N+1),1);
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        if e4s(s4e(j),2)==0
            inddb2((s4e(j)-1)*(N+1)+(1:N+1))=(s4e(j)-1)*(N+1)+(1:N+1);
        end
    end
end
inddb2=nonzeros(inddb2)';

fns = setdiff(1:size(e4s,1)*(N+1), inddb2);

xr = (c4n(n4e(:,2),1)-c4n(n4e(:,1),1))/2;
yr = (c4n(n4e(:,2),2)-c4n(n4e(:,1),2))/2;
xs = (c4n(n4e(:,3),1)-c4n(n4e(:,1),1))/2;
ys = (c4n(n4e(:,3),2)-c4n(n4e(:,1),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

a=a(0);

% d = zeros(size(c4n,1),1);
Fr=zeros((N+1),(N+1),size(n4e,1));
Gr=zeros((N+1),size(n4e,1));
AFr=zeros((N+1)*(N+2)/2,size(n4e,1));
ABr=zeros((N+1)*(N+2)/2,(N+1),size(n4e,1));
V=zeros(size(e4s,1)*(N+1),1);

for j=1:size(n4e,1)
    Ab=zeros((N+1)*(N+2)/2);
    Ac=zeros((N+1)*(N+2)/2);
    Ad=zeros((N+1)*(N+2)/2);

    % Ba=zeros((N+1)*(N+2)/2,(N+1));
    Bb=zeros((N+1)*(N+2)/2,(N+1));

    % Ca=zeros((N+1),(N+1)*(N+2)/2);
    Cb=zeros((N+1),(N+1)*(N+2)/2);

    Aa=e*J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds)+...
              b*J(j)*M2D+...
              -J(j)*((rx(j)*Dr'+sx(j)*Ds')*M2D*a(1)+(ry(j)*Dr'+sy(j)*Ds')*M2D*a(2));
    da=J(j)*M2D*f(c4n(ind4e(j,:),:));

    ht=norm(c4n(n4e(j,3),:)-c4n(n4e(j,2),:));
    n=normal(c4n(n4e(j,3),:)-c4n(n4e(j,2),:));
    h=ht;

    v=a*n;
    
    Da=e*(h/2)*(k/ht)*M;

    Ab(:,en(s4e(j),:,1))=Ab(:,en(s4e(j),:,1))-S*e*(h/2)*((rx(j)*Dr(en(s4e(j),:,1),:)'*M+sx(j)*Ds(en(s4e(j),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j),:,1),:)'*M+sy(j)*Ds(en(s4e(j),:,1),:)'*M)*n(2));
    Ac(en(s4e(j),:,1),:)=Ac(en(s4e(j),:,1),:)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j),:,1),:)+sx(j)*M*Ds(en(s4e(j),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j),:,1),:)+sy(j)*M*Ds(en(s4e(j),:,1),:))*n(2));
    Ad(en(s4e(j),:,1),en(s4e(j),:,1))=Ad(en(s4e(j),:,1),en(s4e(j),:,1))+e*(h/2)*(k/ht)*M;
    
    if e4s(s4e(j),2)~=j
        Ba=S*e*(h/2)*((rx(j)*Dr(en(s4e(j),:,1),:)'*M+sx(j)*Ds(en(s4e(j),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j),:,1),:)'*M+sy(j)*Ds(en(s4e(j),:,1),:)'*M)*n(2));
        Bb(en(s4e(j),:,1),:)=Bb(en(s4e(j),:,1),:)-e*(h/2)*(k/ht)*M;

        Ca=e*(h/2)*((rx(j)*M*Dr(en(s4e(j),:,1),:)+sx(j)*M*Ds(en(s4e(j),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j),:,1),:)+sy(j)*M*Ds(en(s4e(j),:,1),:))*n(2));
        Cb(:,en(s4e(j),:,1))=Cb(:,en(s4e(j),:,1))-e*(h/2)*(k/ht)*M;
    else
        Ba=fliplr(S*e*(h/2)*((rx(j)*Dr(en(s4e(j),:,1),:)'*M+sx(j)*Ds(en(s4e(j),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j),:,1),:)'*M+sy(j)*Ds(en(s4e(j),:,1),:)'*M)*n(2)));
        Bb(en(s4e(j),:,1),:)=Bb(en(s4e(j),:,1),:)-fliplr(e*(h/2)*(k/ht)*M);

        Ca=flip(e*(h/2)*((rx(j)*M*Dr(en(s4e(j),:,1),:)+sx(j)*M*Ds(en(s4e(j),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j),:,1),:)+sy(j)*M*Ds(en(s4e(j),:,1),:))*n(2)));
        Cb(:,en(s4e(j),:,1))=Cb(:,en(s4e(j),:,1))-flip(e*(h/2)*(k/ht)*M);
    end

    if e4s(s4e(j),2)==0
        da=da-S*e*(h/2)*((rx(j)*Dr(en(s4e(j),:,1),:)'*M+sx(j)*Ds(en(s4e(j),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j),:,1),:)'*M+sy(j)*Ds(en(s4e(j),:,1),:)'*M)*n(2))*ue(c4n(ind4s(s4e(j),:,1),:));
        da(en(s4e(j),:,1))=da(en(s4e(j),:,1))-(-e*(h/2)*(k/ht)*M)*ue(c4n(ind4s(s4e(j),:,1),:));
    end


    if v>0
        Ad(en(s4e(j),:,1),en(s4e(j),:,1))=Ad(en(s4e(j),:,1),en(s4e(j),:,1))+(h/2)*(M*v);

        if e4s(s4e(j),2)==j
            Cb(:,en(s4e(j),:,1))=Cb(:,en(s4e(j),:,1))-flip((h/2)*(M*v));
        else
            Cb(:,en(s4e(j),:,1))=Cb(:,en(s4e(j),:,1))-(h/2)*(M*v);
        end
    else
        Da=Da-(h/2)*(M*v);

        if e4s(s4e(j),2)==j
            Bb(en(s4e(j),:,1),:)=Bb(en(s4e(j),:,1),:)+fliplr((h/2)*(M*v));
        else
            Bb(en(s4e(j),:,1),:)=Bb(en(s4e(j),:,1),:)+(h/2)*(M*v);
        end

        if e4s(s4e(j),2)==0
            da(en(s4e(j),:,1))=da(en(s4e(j),:,1))-(h/2)*(M*v)*ue(c4n(ind4s(s4e(j),:,1),:));
        end
    end

    Ar=Aa+Ab+Ac+Ad;
    Br=Ba+Bb;
    Cr=Ca+Cb;

    Fr(:,:,j)=-Cr*(Ar\Br)+Da;
    Gr(:,j)=-Cr*(Ar\da);

    AFr(:,j)=Ar\da;
    ABr(:,:,j)=Ar\Br;

    % d(ind4e(j,:))=d(ind4e(j,:))+da;

end

T=(s4e-1)*(N+1)+(1:N+1);

ind=ind4e';
TA=T';

If=repmat(T,1,(N+1))';
Jf=(repmat(TA(:),1,(N+1)))';
F=sparse(If(:),Jf(:),Fr(:));

Ig=TA;
Jg=ones((N+1),size(s4e,1));
G=sparse(Ig(:),Jg(:),Gr(:));

V(fns)=F(fns,fns)\G(fns);

Iaf=ind;
Jaf=ones((N+1)*(N+2)/2,size(s4e,1));
AF=sparse(Iaf(:),Jaf(:),AFr(:));

Iab=repmat(ind4e,1,(N+1))';
Jab=(repmat(TA(:),1,(N+1)*(N+2)/2))';
AB=sparse(Iab(:),Jab(:),ABr(:));

u=AF-AB*V;

% T = delaunay(c4n(:,1),c4n(:,2));
% trisurf(T,c4n(:,1),c4n(:,2),u,"Linestyle","none")

plot3(c4n(:,1),c4n(:,2),u,'.')