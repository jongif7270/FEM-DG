f = @(x) Simplex2DP(x(:,1),x(:,2),2,3);

N=4;

[x,y]=Nodes2D_equi(10*N);
u=f([x,y]);


[int_x,int_y]=Nodes2D_equi(N);
[int_u]=f([int_x,int_y]);
[r,s]=xytors(int_x,int_y);
[a, b] = rstoab(r, s);
[V]=Vandermonde2D(N,r,s);

num=length(int_x);

W=V\int_u;

basis_ind=zeros((N+1)*(N+2)/2,2);
sk=1;
for i=0:N
    for j=0:N-i
        basis_ind(sk,:)=[i,j];
        sk=sk+1;
    end
end

int_poly=@(x) arrayfun(@(j) Simplex2DP(x(:,1),x(:,2),basis_ind(j,1),basis_ind(j,2)), 1:(N+1)*(N+2)/2)*W;


[ra,sa]=xytors(x,y);
[aa, ba] = rstoab(ra, sa);

z=zeros(1,length(aa));

for i=1:length(aa)
    z(i)=int_poly([aa(i),ba(i)]);
end

plot3(int_x,int_y,int_u,'o')
hold on
plot3(x,y,z,'x')