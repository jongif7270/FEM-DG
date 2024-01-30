function [ind4s] = indexforHDG(xl,xr,yl,yr,Mx,My,k)

[~,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,k);

e4s = computeE4s(n4e);

ind4s = zeros(size(e4s,1),k+1);

for j=1:size(e4s,1)
    ind4s(j,:)=(j-1)*(k+1)+(1:k+1);
end


