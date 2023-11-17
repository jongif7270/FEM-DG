function [v] = rstov(v1,v2,v3)
[x,y] = Nodes2D_equi(N);
[r,s] = xytors(x,y);
v = (r(:)+1)/2*v1+(s(:)+1)/2*v2-(r(:)+s(:))/2*v3;
