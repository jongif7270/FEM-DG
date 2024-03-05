function [c4n,n4e,ind4e,inddb] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,k)
ind4e = zeros(2*Mx*My,(k+1)*(k+2)/2);
tmp = repmat((1:k:k*Mx)',1,My) ...
+ repmat((0:k*(k*Mx+1):((k*Mx+1)*((My-1)*k+1)-1)),Mx,1);
tmp2 = repmat(((k+1):k:(k*Mx+1))',1,My) ...
+ repmat((k*(k*Mx+1):k*(k*Mx+1):((k*Mx+1)*(k*My))),Mx,1);
tmp = tmp(:); tmp2=tmp2(:);
for j=1:k+1
ind4e(1:2:2*Mx*My,(1+(j-1)*(k+2-j/2))+(0:(k+1-j))) = ...
repmat(tmp,1,k+2-j)+repmat(((j-1)*(Mx*k+1)+(0:(k+1-j))),Mx*My,1);
ind4e(2:2:2*Mx*My,(1+(j-1)*(k+2-j/2))+(0:(k+1-j))) = ...
repmat(tmp2,1,k+2-j)+repmat((-(j-1)*(Mx*k+1)-(0:(k+1-j))),Mx*My,1);
end
n4e = ind4e(:,[k+1 (k+1)*(k+2)/2 1]);
inddb = [1:(k*Mx+1), 2*(k*Mx+1):(k*Mx+1):(k*Mx+1)*(k*My+1), ...
((k*Mx+1)*(k*My+1)-1):-1:(k*My*(k*Mx+1)+1), ...
((k*My-1)*(k*Mx+1)+1):-(k*Mx+1):(k*Mx+2)];
%x=linspace(xl,xr,k*Mx+1);
%y=linspace(yl,yr,k*My+1);
%%
[r] = JacobiGL(0,0,k);
hx=(xr-xl)/Mx;
hy=(yr-yl)/My;
x=zeros(1,k*Mx+1);
y=zeros(1,k*My+1);
for j=0:(Mx-1)
    x(j*k+1:j*k+k+1)=xl+j*hx+hx*(r+1)/2;
    y(j*k+1:j*k+k+1)=yl+j*hy+hy*(r+1)/2;
end
%%
y=repmat(y,k*Mx+1,1);
x=repmat(x,k*My+1,1)';
c4n = [x(:), y(:)];
end