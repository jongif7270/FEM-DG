function [ind4e,n4e,inddb,c4n,e4s] = indexforDG2(xl,xr,yl,yr,Mx,My,k)

n4e=zeros(2*Mx*My,3);
c4n=zeros(Mx*My*(k+1)*(k+2),2);
inddb=[];

ind4e=zeros(2*Mx*My,(k+1)*(k+2)/2);
for j=1:2*Mx*My %number of element
    ind4e(j,:)=ind4e(j,:)+(1+(j-1)*(k+1)*(k+2)/2):j*(k+1)*(k+2)/2;
    n4e(j,:)=n4e(j,:)+[ind4e(j,1) ind4e(j,k+1) ind4e(j,(k+1)*(k+2)/2)];
end

%[x,y] = Nodes2D_equi(k);
[x,y] = Nodes2D(k);
[r,s] = xytors(x,y);
[oldc4n,oldn4e,~,~] = mesh_fem_2d_triangle2(xl,xr,yl,yr,Mx,My,k);
for j=1:2*Mx*My
    c4e = (r(:)+1)/2*oldc4n(oldn4e(j,1),:)+(s(:)+1)/2*oldc4n(oldn4e(j,2),:)-(r(:)+s(:))/2*oldc4n(oldn4e(j,3),:);
    jth=ind4e(j,:);
    for i=1:(k+1)*(k+2)/2
        c4n(jth(i),:)=c4e(i,:);
    end
end

a=zeros(1,k+1);
a(1)=k+1;
for j=1:k
    a(j+1)=a(j)+k+1-j;
end

b=zeros(1,k+1);
b(1)=1;
for j=1:k
    b(j+1)=b(j)+k+2-j;
end

e4s = computeE4s(oldn4e);
tmp=[];
for j=2:size(e4s,1)
    if e4s(j,1)<e4s(j-1,1)
        tmp=[tmp j-1];
    end
end

for j=1:tmp(1)
    if e4s(j,2)==0
        inddb=[inddb (a+(e4s(j,1)-1)*(k+1)*(k+2)/2)];
    end
end

for j=tmp(1)+1:tmp(2)
    if e4s(j,2)==0
        inddb=[inddb (b+(e4s(j,1)-1)*(k+1)*(k+2)/2)];
    end
end

for j=tmp(2)+1:size(e4s,1)
    if e4s(j,2)==0
        inddb=[inddb ((1:k+1)+(e4s(j,1)-1)*(k+1)*(k+2)/2)];
    end
end

inddb=unique(inddb);