function[u]=FEM2D[f,N]
u=zeros((N+2)*(N+1)/2,1);
f_value=zeros((N+2)*(N+1)/2,1);
%f_value 찾기
%stiffnessmatrix S 찾기
u=S\M*f_value;