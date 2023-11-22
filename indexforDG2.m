function [ind4e,n4e,c4n] = indexforDG2(xl,xr,yl,yr,Mx,My,k)

n4e=zeros(2*Mx*My,3);
c4n=zeros(Mx*My*(k+1)*(k+2),2);
inddb=[];

ind4e=zeros(2*Mx*My,(k+1)*(k+2)/2);
for j=1:2*Mx*My %number of element
    ind4e(j,:)=ind4e(j,:)+(1+(j-1)*(k+1)*(k+2)/2):j*(k+1)*(k+2)/2;
    n4e(j,:)=n4e(j,:)+[ind4e(j,1) ind4e(j,k+1) ind4e(j,(k+1)*(k+2)/2)];
end

[x,y] = Nodes2D_equi(k);
[r,s] = xytors(x,y);
[oldc4n,oldn4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,k);
for j=1:2*Mx*My
    c4e = (r(:)+1)/2*oldc4n(oldn4e(j,1),:)+(s(:)+1)/2*oldc4n(oldn4e(j,2),:)-(r(:)+s(:))/2*oldc4n(oldn4e(j,3),:);
    jth=ind4e(j,:);
    for i=1:(k+1)*(k+2)/2
        c4n(jth(i),:)=c4e(i,:);
    end
end