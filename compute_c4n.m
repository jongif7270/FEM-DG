function [newind4e,newc4n] = compute_c4n(Mx,My,k)
%%임시%%
[c4n,n4e,~,~] = mesh_fem_2d_triangle(0,1,0,1,Mx,My,1);
%%%%%%%
allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
[~,ind] = unique(sort(allSides,2),'rows','first');
number_of_sides=length(ind);
number_of_elements=length(n4e(:,1));
number_of_nodes=length(c4n(:,1));
ind4e=compute_ind4e(c4n,n4e,k);
newind4e=[ind4e(:,3) ind4e(:,2*k+2:3*k) ind4e(:,1)];
L=3*k+1;
for j=1:k-1
    R=L+k-2-j;
    newind4e=[newind4e ind4e(:,2*k+2-j) ind4e(:,L:R) ind4e(:,3+j)];
    L=R+1;
end
newind4e = [newind4e ind4e(:,2)];
[x,y] = Nodes2D_equi(k);
[r,s] = xytors(x,y);
newc4n=zeros(number_of_nodes+number_of_sides*(k-1)+number_of_elements*((k-2)*(k-1)/2),2);
for j=1:number_of_elements
    c4e = (r(:)+1)/2*c4n(n4e(j,1),:)+(s(:)+1)/2*c4n(n4e(j,2),:)-(r(:)+s(:))/2*c4n(n4e(j,3),:);
    jth=newind4e(j,:);
    for i=1:(k+1)*(k+2)/2
        newc4n(jth(i),:)=c4e(i,:);
    end
end
%%체크%%
trimesh(n4e,c4n(:,1),c4n(:,2))
hold on
for i = 1:size(newc4n, 1)
    text(newc4n(i, 1), newc4n(i, 2), num2str(i));
end

%%%%%%