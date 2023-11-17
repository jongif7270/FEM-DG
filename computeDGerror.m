function err=computeDGerror(c4n,n4e,ind4e,Dr,Ds,u,ux,uy,V,N)
err = 0;
xr = (c4n(n4e(:,1),1)-c4n(n4e(:,3),1))/2;
yr = (c4n(n4e(:,1),2)-c4n(n4e(:,3),2))/2;
xs = (c4n(n4e(:,2),1)-c4n(n4e(:,3),1))/2;
ys = (c4n(n4e(:,2),2)-c4n(n4e(:,3),2))/2;
J = xr.*ys-xs.*yr;
rx=ys./J; ry=-xs./J; sx=-yr./J; sy=xr./J;
I = eye((N+1)*(N+2)/2);
for j=1:size(n4e,1)
M=J(j)*I/(V*V');
Dx_u = (rx(j)*Dr+sx(j)*Ds)*u(ind4e(j,:));
Dy_u = (ry(j)*Dr+sy(j)*Ds)*u(ind4e(j,:));
Dex=ux(c4n(ind4e(j,:),:)) - Dx_u;
Dey=uy(c4n(ind4e(j,:),:)) - Dy_u;
err=err+J(j)*(Dex'*M*Dex+Dey'*M*Dey);
end
err=sqrt(err);
end