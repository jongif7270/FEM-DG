function [u,V2D,Dr,Ds,c4n2] = HDG10(M,N)

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [0.8,0.6];b=0;e=1;S=-1;k=1;  f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +0.8*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) +0.6*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

%% Ex 3.3
% xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [0.8,0.6];b=1;e=0.01;S=1;k=4*N^2;  f=@(x) b.*(1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)) +0.8*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8 +0.6*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4 +e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

% xl=-4;xr=4;yl=-4;yr=4;Mx=M;My=M;    a=@(x) [(x(:,1).^2).*x(:,2), -x(:,1).*(x(:,2).^2)];b=1;e=0.01;S=1;k=4*N^2;  f=@(x) b.*(1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8)) +(x(:,1).^2).*x(:,2).*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,2)+1).^2/8 +-x(:,1).*(x(:,2).^2).*cos(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1).*(x(:,2)+1)/4 +e.*(sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

%% Ex 3.4
% xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))), -exp(x(:,1)).*x(:,2).*sin(x(:,2))]; b=0;e=1;S=-1;k=1; f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +pi*exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2))).*cos(pi*x(:,1)).*sin(pi*x(:,2)) +pi*(-exp(x(:,1)).*x(:,2).*sin(x(:,2))).*sin(pi*x(:,1)).*cos(pi*x(:,2)) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;    a=@(x) [(x(:,1).^2).*x(:,2), -x(:,1).*(x(:,2).^2)]; b=0;e=1;S=-1;k=1; f=@(x) b*sin(pi*x(:,1)).*sin(pi*x(:,2)) +(x(:,1).^2).*x(:,2).*(pi*cos(pi*x(:,1)).*sin(pi*x(:,2))) -x(:,1).*(x(:,2).^2).*(pi*sin(pi*x(:,1)).*cos(pi*x(:,2))) +e*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

int=1;

[c4n,n4e,~,~] = mesh_fem_2d_triangle2(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,~,c4n2,~] = indexforDG2(xl,xr,yl,yr,Mx,My,N);
[~,in4e,~,~] = mesh_fem_2d_triangle2(xl,xr,yl,yr,Mx,My,int*N);
[iind4e,~,~,ic4n2,~] = indexforDG2(xl,xr,yl,yr,Mx,My,int*N);

[r1D] = JacobiGL(0,0,N);
[int_r1D] = JacobiGL(0,0,int*N);
[V1D] = Vandermonde1D(N,r1D);
[int_V1D]=Vandermonde1D(N,int_r1D);
[iV1D]=Vandermonde1D(int*N,int_r1D);
I1D = eye(size(V1D,1));
iI1D = eye(size(iV1D,1));

[x,y]=Nodes2D(N);
[int_x,int_y]=Nodes2D(int*N);
[r,s]=xytors(x,y);
[int_r,int_s]=xytors(int_x,int_y);
V2D=Vandermonde2D(N,r,s);
int_V2D=Vandermonde2D(N,int_r,int_s);
iV2D=Vandermonde2D(int*N,int_r,int_s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
[iDr,iDs]=Dmatrices2D(int*N,int_r,int_s,iV2D);
I2D = eye((N+1)*(N+2)/2);
iI2D = eye((int*N+1)*(int*N+2)/2);

A1D=V1D'\int_V1D';
A2D=V2D'\int_V2D';

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

iM=iI1D/(iV1D*iV1D');
iM2D=iI2D/(iV2D*iV2D');

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
s4e = computeS4e(n4e);
en=mod(ind4s(:,:,:)-1,(N+1)*(N+2)/2)+1;

iind4s=edge(in4e,int*N);
ie4s = computeE4s(in4e);
is4e = computeS4e(in4e);

% inddb1=zeros(size(c4n2,1),1);
% for i=1:size(ind4s,1)
%     if ind4s(i,:,1)==ind4s(i,:,2)
%         inddb1(ind4s(i,:,1))=ind4s(i,:,1);
%     end
% end
% inddb1=nonzeros(sort(inddb1))';
% 
% fns2 = setdiff(1:size(c4n2,1),inddb1);

inddb2=zeros(size(e4s,1)*(N+1),1);
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        if e4s(s4e(j,i),2)==0
            inddb2((s4e(j,i)-1)*(N+1)+(1:N+1))=(s4e(j,i)-1)*(N+1)+(1:N+1);
        end
    end
end
inddb2=nonzeros(inddb2)';

fns = setdiff(1:size(e4s,1)*(N+1), inddb2);

xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

ae=zeros(size(iind4e,2),2,size(iind4e,1));
for j=1:size(iind4e,1)
    for l=1:(int*N+1)*(int*N+2)/2
        ae(l,:,j)=a(ic4n2(((int*N+1)*(int*N+2)/2)*(j-1)+l,:));
    end
    % ae(:,:,j)=a(ic4n2(iind4e(j,:),:));
end

al=zeros(int*N+1,2,size(ie4s,1));
for j=1:size(in4e,1)
    for i=1:size(in4e,2)
        for l=1:int*N+1
            al(l,:,is4e(j,i))=a(ic4n2(iind4s(is4e(j,i),l,1),:));
        end
        % al(:,:,is4e(j,i))=a(ic4n2(iind4s(is4e(j,i),:,1),:));
    end
end


% d = zeros(size(ind4e(:),1),1);
Da=zeros(3*(N+1),3*(N+1));
Fr=zeros(3*(N+1),3*(N+1),size(n4e,1));
Gr=zeros(3*(N+1),size(n4e,1));
AFr=zeros((N+1)*(N+2)/2,size(n4e,1));
ABr=zeros((N+1)*(N+2)/2,3*(N+1),size(n4e,1));
V=zeros(size(e4s,1)*(N+1),1);

for j=1:size(n4e,1)
    Ab=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);
    Ac=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);
    Ad=zeros((N+1)*(N+2)/2,(N+1)*(N+2)/2);

    Ba=zeros((N+1)*(N+2)/2,3*(N+1));
    Bb=zeros((N+1)*(N+2)/2,3*(N+1));

    Ca=zeros(3*(N+1),(N+1)*(N+2)/2);
    Cb=zeros(3*(N+1),(N+1)*(N+2)/2);

    dx=rx(j)*iDr'+sx(j)*iDs';
    dy=ry(j)*iDr'+sy(j)*iDs';
    a1=ae(:,1,j);
    a2=ae(:,2,j);

    % a1=a1'; a2=a2';

    Aa=e*J(j)*((rx(j)^2+ry(j)^2)*Dr'*M2D*Dr+(rx(j)*sx(j)+ry(j)*sy(j))*(Ds'*M2D*Dr+Dr'*M2D*Ds)+(sx(j)^2+sy(j)^2)*Ds'*M2D*Ds)+...
              b*J(j)*M2D+...
              -J(j)*A2D*(a1.*(dx*iM2D)+a2.*(dy*iM2D))*A2D';
    da=J(j)*M2D*f(c4n2(ind4e(j,:),:));

    ht=norm(c4n(n4e(j,2),:)-c4n(n4e(j,1),:));
    for i=1:3
        n=normal(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));
        h=norm(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));

        v=al(:,:,is4e(j,i))*n;

        % v=v';

        % vp=0;
        % vn=0;
        % for l=1:int*N+1
        %     if v(l)>=0
        %         vp=vp+1;
        %     elseif v(l)<0
        %         vn=vn+1;
        %     end
        % end

        % vm=[a1(ic4n2(iind4s(s4e(j,i),ceil((N+1)/2),1),:)) a2(ic4n2(iind4s(s4e(j,i),ceil((N+1)/2),1),:))]*n;
        % v1=[a1(c4n2(ind4s(s4e(j,i),:,1),:)) a2(c4n2(ind4s(s4e(j,i),:,l),:))]*n;
        % v2=[a1(c4n2(ind4s(s4e(j,i),:,2),:)) a2(c4n2(ind4s(s4e(j,i),:,2),:))]*n;

        if v(ceil((N*int+1)/2))>0
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=e*(h/2)*(k/ht)*M;
    
            Ab(:,en(s4e(j,i),:,1))=Ab(:,en(s4e(j,i),:,1))-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:)=Ac(en(s4e(j,i),:,1),:)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
    
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))+e*(h/2)*(k/ht)*M+(h/2)*A1D*(iM.*v')*A1D';
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M-(h/2)*A1D*(iM.*v')*A1D';
            else
                Ba(:,(i-1)*(N+1)+(1:N+1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))-e*(h/2)*(k/ht)*M;
    
                Ca((i-1)*(N+1)+(1:N+1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M-(h/2)*A1D*(iM.*v')*A1D';
            end
            if e4s(s4e(j,i),2)==0
                da=da-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2))*ue(c4n2(ind4s(s4e(j,i),:,1),:));
                da(en(s4e(j,i),:,1))=da(en(s4e(j,i),:,1))-(-e*(h/2)*(k/ht)*M)*ue(c4n2(ind4s(s4e(j,i),:,1),:));
            end
        else
            Da((i-1)*(N+1)+(1:N+1),(i-1)*(N+1)+(1:N+1))=e*(h/2)*(k/ht)*M-(h/2)*A1D*(iM.*v')*A1D';
    
            Ab(:,en(s4e(j,i),:,1))=Ab(:,en(s4e(j,i),:,1))-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
            Ac(en(s4e(j,i),:,1),:)=Ac(en(s4e(j,i),:,1),:)-e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
            Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))=Ad(en(s4e(j,i),:,1),en(s4e(j,i),:,1))+e*(h/2)*(k/ht)*M;
    
            if e4s(s4e(j,i),2)==j
                Ba(:,(i-1)*(N+1)+(N+1:-1:1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(N+1:-1:1))-e*(h/2)*(k/ht)*M+(h/2)*A1D*(iM.*flip(v)')*A1D';
    
                Ca((i-1)*(N+1)+(N+1:-1:1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(N+1:-1:1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M;
            else
                Ba(:,(i-1)*(N+1)+(1:N+1))=S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2));
                Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))=Bb(en(s4e(j,i),:,1),(i-1)*(N+1)+(1:N+1))-e*(h/2)*(k/ht)*M+(h/2)*A1D*(iM.*v')*A1D';
    
                Ca((i-1)*(N+1)+(1:N+1),:)=e*(h/2)*((rx(j)*M*Dr(en(s4e(j,i),:,1),:)+sx(j)*M*Ds(en(s4e(j,i),:,1),:))*n(1)+(ry(j)*M*Dr(en(s4e(j,i),:,1),:)+sy(j)*M*Ds(en(s4e(j,i),:,1),:))*n(2));
                Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))=Cb((i-1)*(N+1)+(1:N+1),en(s4e(j,i),:,1))-e*(h/2)*(k/ht)*M;
            end

            if e4s(s4e(j,i),2)==0
                da=da-S*e*(h/2)*((rx(j)*Dr(en(s4e(j,i),:,1),:)'*M+sx(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(1)+(ry(j)*Dr(en(s4e(j,i),:,1),:)'*M+sy(j)*Ds(en(s4e(j,i),:,1),:)'*M)*n(2))*ue(c4n2(ind4s(s4e(j,i),:,1),:));
                da(en(s4e(j,i),:,1))=da(en(s4e(j,i),:,1))-(-e*(h/2)*(k/ht)*M+(h/2)*A1D*(iM.*v')*A1D')*ue(c4n2(ind4s(s4e(j,i),:,1),:));
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

    % d(ind4e(j,:))=d(ind4e(j,:))+da;
end

T=[(s4e(:,1)-1)*(N+1)+(1:N+1) (s4e(:,2)-1)*(N+1)+(1:N+1) (s4e(:,3)-1)*(N+1)+(1:N+1)];

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

% u=zeros(size(c4n2,1),1);
% u(fns2)=AF(fns2)-AB(fns2,:)*V;

% plot3(c4n2(:,1),c4n2(:,2),u,'.')
T = delaunay(c4n2(:,1),c4n2(:,2));
trisurf(T,c4n2(:,1),c4n2(:,2),u,"Linestyle","none")