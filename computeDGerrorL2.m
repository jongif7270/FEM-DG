function err=computeDGerrorL2(c4n,c4n2,n4e,ind4e,Dr,Ds,u,ux,uy,V,N)
err = 0;
xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
I = eye((N+1)*(N+2)/2);
M=I/(V*V');

uo=@(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

for j=1:size(n4e,1)
ue=uo(c4n2(ind4e(j,:),:)) - u(ind4e(j,:));
err=err+J(j)*(ue'*M*ue);
end
err=sqrt(err);
end