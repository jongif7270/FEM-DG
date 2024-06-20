function [u,V2D,Dr,Ds,c4n] = FEMHDG3(M,N)


% xl=-4;xr=4;yl=-4;yr=4;Mx=M;My=M;

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

% a=@(x) [0,0];b=0;e=1;S=1;k=4*N^2;  f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +0*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) +0*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% a=@(x) [0.8,0.6];b=1;e=1;S=1;k=4*N^2;  f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +0.8*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) +0.6*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% a=@(x) [0,0];b=0;e=1;S=1;k=4*N^2;  f=@(x) b.*(1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)) +0*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8 +0*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4 +e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

% a=@(x) [0.8,0.6];b=1;e=1;S=1;k=4*N^2;  f=@(x) b.*(1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)) +0.8*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8 +0.6*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4 +e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

a=@(x) [exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))), -exp(x(:,1)).*x(:,2).*sin(x(:,2))]; b=0;e=1;S=-1;k=1; f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +pi*exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))).*cos(pi*x(:,1)).*sin(pi*x(:,2)) +pi*(-exp(x(:,1)).*x(:,2).*sin(x(:,2))).*sin(pi*x(:,1)).*cos(pi*x(:,2)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

[~,n4e,ind4e,~,ind4s,e4s,s4e,~,en,ind4p,s4p,e4p,c4n] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

%%
int=3;
[~,n4e2,ind4e2,~,~,~,~,~,~,~,~,~,c4n2] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N*int);

[r1D2] = JacobiGL(0,0,N*int);
% [r1D] = Nodes1D_equi(N);
[V1D2] = Vandermonde1D(N*int,r1D2);
[iV1D] = Vandermonde1D(N,r1D2);
I1D2 = eye(size(V1D2,1));

[x2,y2]=Nodes2D(N*int);
% [x,y]=Nodes2D_equi(N);
[r2,s2]=xytors(x2,y2);
V2D2=Vandermonde2D(N*int,r2,s2);
iV2D=Vandermonde2D(N,r2,s2);
[Dr2,Ds2]=Dmatrices2D(N*int,r2,s2,V2D2);
I2D2 = eye((N*int+1)*(N*int+2)/2);

M2=I1D2/(V1D2*V1D2');
M2D2=I2D2/(V2D2*V2D2');


ae = zeros(size(ind4e2,2),2,size(ind4e2,1));
for j=1:size(ind4e2,1)
    for l=1:size(ind4e2,2)
        ae(l,:,j)=a(c4n2(ind4e2(j,l),:));
    end
end

al=zeros(int*N+1,2,size(n4e2,1));
for j=1:size(n4e2,1)
    for i=1:int*N+1
        al(i,:,j)=a(c4n2(ind4s(s4e(j),:,1)));
    end
end

%%

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

A1D=V1D'\iV1D';
A2D=V2D'\iV2D';

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

Fr=zeros(4*(N+1),4*(N+1),size(ind4p,1));
Gr=zeros(4*(N+1),size(ind4p,1));
AFr=zeros(size(ind4p,2),size(ind4p,1));
ABr=zeros(size(ind4p,2),4*(N+1),size(ind4p,1));
V=zeros(size(e4s,1)*(N+1),1);

for j = 1:size(ind4p,1)
    Aa=zeros(size(ind4p,2));
    Ab=Aa;
    Ac=Aa;
    Ad=Aa;

    Ba=zeros(size(ind4p,2),4*(N+1));
    Bb=Ba;

    Ca=zeros(4*(N+1),size(ind4p,2));
    Cb=Ca;

    Da=zeros(4*(N+1));

    da=zeros(size(ind4p,2),1);
    for i = 1:4
        ht=norm(c4n(n4e(e4p(j,i),3),:)-c4n(n4e(e4p(j,i),2),:));
        h=ht;
        n=normal(c4n(n4e(e4p(j,i),3),:)-c4n(n4e(e4p(j,i),2),:));
        v=al(:,:,e4p(j,i))*n;

        da(ind4e(e4p(1,i),:)) = da(ind4e(e4p(1,i),:)) + J(j)*M2D*f(c4n(ind4e(e4p(j,i),:),:));

        vu = (h/2)*M;
        vu2 = (h/2)*M2;
        gvu = (h/2)*((rx(e4p(j,i))*Dr(en(s4e(e4p(j,i)),:,1),:)'*M+sx(e4p(j,i))*Ds(en(s4e(e4p(j,i)),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(e4p(j,i)),:,1),:)'*M+sy(e4p(j,i))*Ds(en(s4e(e4p(j,i)),:,1),:)'*M)*n(2));
        vgu = (h/2)*((rx(e4p(j,i))*M*Dr(en(s4e(e4p(j,i)),:,1),:)+sx(e4p(j,i))*M*Ds(en(s4e(e4p(j,i)),:,1),:))*n(1)+(ry(e4p(j,i))*M*Dr(en(s4e(e4p(j,i)),:,1),:)+sy(e4p(j,i))*M*Ds(en(s4e(e4p(j,i)),:,1),:))*n(2));


        Aa(ind4e(i,:),ind4e(i,:))=Aa(ind4e(i,:),ind4e(i,:))+...
            e*J(e4p(j,i))*((rx(e4p(j,i))^2+ry(e4p(j,i))^2)*Dr'*M2D*Dr+(rx(e4p(j,i))*sx(e4p(j,i))+ry(e4p(j,i))*sy(e4p(j,i)))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(e4p(j,i))^2+sy(e4p(j,i))^2)*Ds'*M2D*Ds)+...
            b*J(e4p(j,i))*M2D+...
            -J(e4p(j,i))*A2D*(ae(:,1,e4p(j,i)).*(rx(e4p(j,i))*Dr2'+sx(e4p(j,i))*Ds2')*M2D2+ae(:,2,e4p(j,i)).*(ry(e4p(j,i))*Dr2'+sy(e4p(j,i))*Ds2')*M2D2)*A2D';

        Ab(ind4e(i,:),ind4s(i,:,1))=Ab(ind4e(i,:),ind4s(i,:,1)) - e*S*gvu;
        Ac(ind4s(i,:,1),ind4e(i,:))=Ac(ind4s(i,:,1),ind4e(i,:)) - e*vgu;
        Ad(ind4s(i,:,1),ind4s(i,:,1))=Ad(ind4s(i,:,1),ind4s(i,:,1)) + e*(k/ht)*vu;

        if e4s(s4e(e4p(j,i)),2)==j
            Ba(ind4e(i,:),(i-1)*(N+1)+(N+1:-1:1))=Ba(ind4e(i,:),(i-1)*(N+1)+(N+1:-1:1)) + e*S*gvu;
            Bb(ind4s(i,:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(ind4s(i,:,1),(i-1)*(N+1)+(N+1:-1:1)) -e*(k/ht)*vu;

            Ca((i-1)*(N+1)+(N+1:-1:1),ind4e(i,:))=Ca((i-1)*(N+1)+(N+1:-1:1),ind4e(i,:)) + e*S*vgu;
            Cb((i-1)*(N+1)+(N+1:-1:1),ind4s(i,:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),ind4s(i,:,1)) -e*(k/ht)*vu;
        else
            Ba(ind4e(i,:),(i-1)*(N+1)+(1:N+1))=Ba(ind4e(i,:),(i-1)*(N+1)+(1:N+1)) + e*S*gvu;
            Bb(ind4s(i,:,1),(i-1)*(N+1)+(1:N+1))=Bb(ind4s(i,:,1),(i-1)*(N+1)+(1:N+1)) -e*(k/ht)*vu;

            Ca((i-1)*(N+1)+(1:N+1),ind4e(i,:))=Ca((i-1)*(N+1)+(1:N+1),ind4e(i,:)) + e*S*vgu;
            Cb((i-1)*(N+1)+(1:N+1),ind4s(i,:,1))=Cb((i-1)*(N+1)+(1:N+1),ind4s(i,:,1)) -e*(k/ht)*vu;
        end

        Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1)) + e*(k/ht)*vu;

        if v(ceil((N*int+1)/2))>0
            Ad(ind4s(i,:,1),ind4s(i,:,1))=Ad(ind4s(i,:,1),ind4s(i,:,1)) + A1D*(vu2.*v')*A1D';

            if e4s(s4e(e4p(j,i)),2)==j
                Cb((i-1)*(N+1)+(N+1:-1:1),ind4s(i,:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),ind4s(i,:,1)) - A1D*(vu2.*v')*A1D';
            else
                Cb((i-1)*(N+1)+(1:N+1),ind4s(i,:,1))=Cb((i-1)*(N+1)+(1:N+1),ind4s(i,:,1)) - A1D*(vu2.*v')*A1D';
            end

            if e4s(s4e(e4p(j,i)),2)==0
                da(ind4e(i,:)) = da(ind4e(i,:))-S*e*gvu*ue(c4n(ind4s(s4e(e4p(j,i)),:,1),:));
                da(ind4s(s4e(i),:,1))=da(ind4s(s4e(i),:,1)) + e*(k/ht)*vu*ue(c4n(ind4s(s4e(e4p(j,i)),:,1),:));
            end
        else
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1)) - A1D*(vu2.*v')*A1D';

            if e4s(s4e(e4p(j,i)),2)==j
                Bb(ind4s(i,:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(ind4s(i,:,1),(i-1)*(N+1)+(N+1:-1:1)) + A1D*(vu2.*flip(v)')*A1D';
            else
                Bb(ind4s(i,:,1),(i-1)*(N+1)+(1:N+1))=Bb(ind4s(i,:,1),(i-1)*(N+1)+(1:N+1)) + A1D*(vu2.*v')*A1D';
            end

            if e4s(s4e(e4p(j,i)),2)==0
                da(ind4e(e4p(1,i),:)) = da(ind4e(e4p(1,i),:))-S*e*gvu*ue(c4n(ind4s(s4e(e4p(j,i)),:,1),:));
                da(ind4s(s4e(e4p(1,i)),:,1))=da(ind4s(s4e(e4p(1,i)),:,1)) + (e*(k/ht)*vu-A1D*(vu2.*v')*A1D')*ue(c4n(ind4s(s4e(e4p(j,i)),:,1),:));
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
end


% T=[(s4p(:,1)-1)*(N+1)+(1:N+1) (s4p(:,2)-1)*(N+1)+(1:N+1) (s4p(:,3)-1)*(N+1)+(1:N+1) (s4p(:,4)-1)*(N+1)+(1:N+1)];

T = zeros(size(s4p,1),(N+1)*4);
for j=1:size(s4p,1)
    for i=1:4
        if e4s(s4p(j,i),2)==e4p(j,i)
            T(j,(i-1)*(N+1)+(1:N+1))=(s4p(j,i)-1)*(N+1)+flip(1:N+1);
        else
            T(j,(i-1)*(N+1)+(1:N+1))=(s4p(j,i)-1)*(N+1)+(1:N+1);
        end
    end
end

ind=ind4p';
TA=T';

If=repmat(T,1,4*(N+1))';
Jf=(repmat(TA(:),1,4*(N+1)))';
F=sparse(If(:),Jf(:),Fr(:));

Ig=TA;
Jg=ones((N+1),size(s4e,1));
G=sparse(Ig(:),Jg(:),Gr(:));

V(fns)=F(fns,fns)\G(fns);

Iaf=ind;
Jaf=ones(size(ind4p,1),size(ind4p,2));
AF=sparse(Iaf(:),Jaf(:),AFr(:));

Iab=repmat(ind4p,1,4*(N+1))';
Jab=(repmat(TA(:),1,size(ind4p,2)))';
AB=sparse(Iab(:),Jab(:),ABr(:));

u=AF-AB*V;

% plot3(c4n(:,1),c4n(:,2),u,'.')
T = delaunay(c4n(:,1),c4n(:,2));
trisurf(T,c4n(:,1),c4n(:,2),u,"Linestyle","none")