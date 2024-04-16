function [c4n,n4e,ind4e,inddb,ind4snew,e4s,s4e,n4s,en,ind4p,s4p,e4p] = mesh_FEMDG(xl,xr,yl,yr,Mx,My,N)
% 
% xl=0;
% xr=1;
% yl=0;
% yr=1;
% Mx=2;
% My=2;
% N=2;

% c4n=zeros(((2*N+1)*(N+1)-N)*Mx*My,2);

% xl=0;xr=1;yl=0;yr=1;Mx=M;My=M;

a=linspace(xl,xr,Mx+1);
b=linspace(yl,yr,My+1);
c4n=[];
for j=1:My
    for i=1:Mx
        x=linspace(a(i),a(i+1),2*N+1);
        y=linspace(b(j),b(j+1),2*N+1);
        x=repmat(x,2*N+1,1)';
        y=repmat(y,2*N+1,1);
        A = [x(:), y(:)];
        even_rows_index=2:2:size(A, 1);
        A= A(setdiff(1:size(A, 1), even_rows_index), :);
        c4n=[c4n;A];
    end
end
% plot(c4n(:,1),c4n(:,2),'.')
% for i = 1:size(c4n, 1)
%     text(c4n(i,1), c4n(i,2), num2str(i));
% end

% ind4e= zeros(4,(N+1)*(N+2)/2,Mx*My);
% 
% for k=1:Mx*My
%     sk=1;
%     for i=0:N
%         for j=0:N-i
%             ind4e(1,sk,k)=N*(N+1)+1 +j*N -i*(N+1) +(k-1)*((2*N+1)*(N+1)-N);
%             ind4e(2,sk,k)=N*(N+1)+1 -j*(N+1) -i*N +(k-1)*((2*N+1)*(N+1)-N);
%             ind4e(3,sk,k)=N*(N+1)+1 -j*N +i*(N+1) +(k-1)*((2*N+1)*(N+1)-N);
%             ind4e(4,sk,k)=N*(N+1)+1 +j*(N+1) +i*N +(k-1)*((2*N+1)*(N+1)-N);
%             sk=sk+1;
%         end
%     end
% end

ind4e= zeros(Mx*My*4,(N+1)*(N+2)/2);

for k=1:Mx*My
    sk=1;
    for i=0:N
        for j=0:N-i
            ind4e(4*(k-1)+1,sk)=N*(N+1)+1 +j*N -i*(N+1) +(k-1)*((2*N+1)*(N+1)-N);
            ind4e(4*(k-1)+2,sk)=N*(N+1)+1 -j*(N+1) -i*N +(k-1)*((2*N+1)*(N+1)-N);
            ind4e(4*(k-1)+3,sk)=N*(N+1)+1 -j*N +i*(N+1) +(k-1)*((2*N+1)*(N+1)-N);
            ind4e(4*(k-1)+4,sk)=N*(N+1)+1 +j*(N+1) +i*N +(k-1)*((2*N+1)*(N+1)-N);
            sk=sk+1;
        end
    end
end


n4e = [ind4e(:,1) ind4e(:,N+1) ind4e(:,(N+1)*(N+2)/2)];

ind4s = zeros(Mx*My*4,N+1);
for k=1:Mx*My
    for j=1:4
        tmp=0;
        for i=1:N+1
            ind4s((k-1)*4+j,i)=ind4e((k-1)*4+j,N+1+tmp);
            tmp=tmp+(N+1-i);
        end
    end
end

tmp=[];
for i=1:Mx
    tmp=[tmp 4*(i-1)+2 4*(My-1)*Mx+4*i];
end
for j=1:My
    tmp=[tmp 4*((j-1)*Mx)+1 4*(j*Mx-1)+3];
end

inddb1=sort(tmp);

a=ind4s(tmp,:,1);
inddb=unique(a(:)');

e4s = zeros(Mx*My*4,2);
for j=1:Mx*My*4
    e4s(j,1)=j;
    if ismember(j,inddb1)==1
        e4s(j,2)=0;
    else
        if mod(j,4)==1
            e4s(j,2)=j-2;
        elseif mod(j,4)==2
            e4s(j,2)=j-4*Mx+2;
        elseif mod(j,4)==3
            e4s(j,2)=j+2;
        else
            e4s(j,2)=j+4*Mx-2;
        end
    end
end

tmp=zeros(Mx*My*4,1);
for j=1:Mx*My*4
    if max(e4s(j,1)<e4s(j,2),e4s(j,2)==0)==1
        tmp(j)=j;
    else
        tmp(j)=0;
    end
end
e4s=e4s(nonzeros(tmp),:);

ind4snew=zeros(size(e4s,1),N+1,2);
for j=1:size(e4s,1)
    if e4s(j,2)==0
        ind4snew(j,:,2)=ind4s(e4s(j,1),:);
    else
        ind4snew(j,:,2)=flip(ind4s(e4s(j,2),:));
    end
end
ind4snew(:,:,1)=ind4s(e4s(:,1),:);

n4s=[ind4snew(:,1,1) ind4snew(:,end,1)];

e=zeros(1,N+1);
tmp=0;
for j=1:N+1
   e(j)=N+1+tmp;
   tmp=tmp+(N+1-j);
end

en=zeros(size(e4s,1),N+1,2);
for j=1:size(e4s,1)
    if e4s(j,2)==0
        en(j,:,2)=e;
    else
        en(j,:,2)=flip(e);
    end
end
en(:,:,1)=repmat(e,size(e4s,1),1);

s4e=zeros(size(n4e,1),1);

for j = 1:size(n4e,1)
    for i = 1:size(e4s,1)
        if isequal(e4s(i,1),j) | isequal(e4s(i,2),j)
            s4e(j)=i;
        end
    end
end

% c4n2=zeros(4*Mx*My*(k+1)*(k+2),2);
% [x,y] = Nodes2D(k);
% [r,s] = xytors(x,y);
% for j=1:2*Mx*My
%     c4e = (r(:)+1)/2*c4n(n4e(j,1),:)+(s(:)+1)/2*c4n(n4e(j,2),:)-(r(:)+s(:))/2*c4n(n4e(j,3),:);
%     jth=ind4e(j,:);
%     for i=1:(k+1)*(k+2)/2
%         c4n2(jth(i),:)=c4e(i,:);
%     end
% end

ind4p = zeros(Mx*My,size(c4n,1)/(Mx*My));
for j = 1:Mx*My
    ind4p(j,:)=(j-1)*size(c4n,1)/(Mx*My) + (1:size(c4n,1)/(Mx*My));
end

s4p = zeros(Mx*My,4);
for j = 1:Mx*My
    s4p(j,:)=s4e((j-1)*4+(1:4));
end

e4p = zeros(Mx*My,4);
for j = 1:Mx*My
    e4p(j,:)=(j-1)*4+(1:4);
end