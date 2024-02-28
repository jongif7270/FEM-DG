f = @(x) sin(2*x)+cos(x);

N=8;

t=linspace(-1,1,5*N+1);
y=f(t);


[int_t] = Nodes1D_equi(N);
[int_y]=f(int_t);
[V] = Vandermonde1D(N,int_t);

num=length(int_t);

W=V\int_y';

int_poly=@(x) arrayfun(@(j) JacobiP(x, 0, 0, j), 0:N)*W;

z=zeros(1,length(t));

for i=1:length(t)
    z(i)=int_poly(t(i));
end

plot(int_t,int_y,'o')
hold on
plot(t,z,'x')