function [u,V2D,Dr,Ds,c4n] = FEMHDG4(M,N)

xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

% S=1;k=4*N^2;  f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)); ue=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

S=1;k=4*N^2;  f=@(x) (sin(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*(pi^2/16*(x(:,2)+1).^2.*((x(:,1)+1).^2+(x(:,2)+1).^2/4))-cos(pi*(x(:,1)+1).*(x(:,2)+1).^2/8).*pi.*(x(:,1)+1)/4); ue=@(x) 1+sin(pi.*(x(:,1)+1).*(x(:,2)+1).^2/8);

[~,n4e,ind4e,~,ind4s,e4s,s4e,~,en,ind4p,s4p,e4p,c4n] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N);

[r1D] = JacobiGL(0,0,N);
[V1D] = Vandermonde1D(N,r1D);
I1D = eye(size(V1D,1));

[x,y]=Nodes2D(N);
[r,s]=xytors(x,y);
V2D=Vandermonde2D(N,r,s);
[Dr,Ds]=Dmatrices2D(N,r,s,V2D);
I2D = eye((N+1)*(N+2)/2);

M=I1D/(V1D*V1D');
M2D=I2D/(V2D*V2D');

xr = (c4n(n4e(:,2),1)-c4n(n4e(:,1),1))/2;
yr = (c4n(n4e(:,2),2)-c4n(n4e(:,1),2))/2;
xs = (c4n(n4e(:,3),1)-c4n(n4e(:,1),1))/2;
ys = (c4n(n4e(:,3),2)-c4n(n4e(:,1),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;

h = vecnorm(c4n(n4e(:,3),:) - c4n(n4e(:,2),:),2,2);

sign = zeros(size(s4e));
for i = 1:size(e4p,1)
    for j = 1:size(e4p,2)
        if e4s(s4p(i,j),2)==0
            sign(e4p(i,j))=1;
        else
            if e4s(s4p(i,j),1)==e4p(i,j)
                sign(e4p(i,j))=1;
            else
                sign(e4p(i,j))=-1;
            end
        end
    end
end

n = ((c4n(n4e(:,3),:) - c4n(n4e(:,2),:))./h)*[0 -1 ; 1 0];
n = sign.*n;

inddb=zeros(size(e4p,1),4*(N+1));
for i = 1:size(e4p,1)
    for j = 1:size(e4p,2)
        if e4s(s4p(i,j),2)==0
            tmp=mod(e4s(s4p(i,j),1)-1,4);
            inddb(i,tmp*(N+1)+(1:N+1))=tmp*(N+1)+(1:N+1);
        end
    end
end

[GUGV,GUV,UGV,UV] = compute_matrix(rx,ry,sx,sy,J,h,n,Dr,Ds,M2D,M,en,sign);

U=zeros(size(c4n,1),1);
for i = 1:size(e4p,1)
    A=zeros(size(ind4p,2),size(ind4p,2));
    B=zeros(size(ind4p,2),4*(N+1));
    C=B';
    D=zeros(4*(N+1),4*(N+1));
    F=zeros(size(ind4p,2),1);
    L=zeros(4*(N+1),1);
    for j = 1:size(e4p,2)
        A(ind4e(j,:),ind4e(j,:))=A(ind4e(j,:),ind4e(j,:))+GUGV(:,:,e4p(i,j));
        A(ind4e(j,:),ind4s(j,:,1))=A(ind4e(j,:),ind4s(j,:,1))-UGV(:,:,e4p(i,j));
        A(ind4s(j,:,1),ind4e(j,:))=A(ind4s(j,:,1),ind4e(j,:))-GUV(:,:,e4p(i,j));
        A(ind4s(j,:,1),ind4s(j,:,1))=A(ind4s(j,:,1),ind4s(j,:,1))+k/h(e4p(i,j))*UV(:,:,e4p(i,j));

        B(ind4e(j,:),(j-1)*(N+1)+(1:N+1))=B(ind4e(j,:),(j-1)*(N+1)+(1:N+1))+S*UGV(:,:,e4p(i,j));
        B(ind4s(j,:,1),(j-1)*(N+1)+(1:N+1))=B(ind4s(j,:,1),(j-1)*(N+1)+(1:N+1))-k/h(e4p(i,j))*UV(:,:,e4p(i,j));

        C((j-1)*(N+1)+(1:N+1),ind4e(j,:))=C((j-1)*(N+1)+(1:N+1),ind4e(j,:))+S*GUV(:,:,e4p(i,j));
        C((j-1)*(N+1)+(1:N+1),ind4s(j,:,1))=C((j-1)*(N+1)+(1:N+1),ind4s(j,:,1))-k/h(e4p(i,j))*UV(:,:,e4p(i,j));

        D((j-1)*(N+1)+(1:N+1),(j-1)*(N+1)+(1:N+1))=D((j-1)*(N+1)+(1:N+1),(j-1)*(N+1)+(1:N+1))+k/h(e4p(i,j))*UV(:,:,e4p(i,j));

        F(ind4e(j,:))=F(ind4e(j,:))+J(e4p(i,j))*M2D*f(c4n(ind4e(e4p(i,j),:),:)); % << dirichlet boundary 설정 / B 와 구성 같음 / ue or g 사용
        if e4s(s4p(i,j),2)==0
            F(ind4e(j,:))=F(ind4e(j,:))-S*UGV(:,:,e4p(i,j))*ue(c4n(ind4s(s4e(e4p(i,j)),:,1),:));
            F(ind4s(j,:,1))=F(ind4s(j,:,1))+k/h(e4p(i,j))*UV(:,:,e4p(i,j))*ue(c4n(ind4s(s4e(e4p(i,j)),:,1),:));
        end
    end
    invA_F = A\F;
    invA_B = A\B;
    fns=setdiff((1:4*(N+1)),nonzeros(inddb(i,:))');
    L1=-(-C*invA_B+D);
    L2=(C*invA_F);
    L(fns) = L1(fns,fns)\L2(fns);
    U(ind4p(i,:))=U(ind4p(i,:))+invA_F-invA_B*L;
end

plot3(c4n(:,1),c4n(:,2),U,'.')

% T = delaunay(c4n(:,1),c4n(:,2));
% trisurf(T,c4n(:,1),c4n(:,2),U,"Linestyle","none")    