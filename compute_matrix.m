function [GUGV,GUV,UGV,UV] = compute_matrix(rx,ry,sx,sy,J,h,n,Dr,Ds,M2D,M)

GUGV=zeros(size(M2D,1),size(M2D,2),size(rx,1));
for i=1:size(rx,1)
    GUGV(:,:,i)=J(i)*(rx(i)^2+ry(i)^2)*Dr'*M2D*Dr+(rx(i)*sx(i)+ry(i)*sy(i))*(Dr'*M2D*Ds+Ds'*M2D*Dr)+(sx(i)^2+sy(i)^2)*Ds'*M2D*Ds;
end

GUV=zeros(size(M,1),size(M2D,2),size(rx,1));
for i=1:size(rx,1)
    GUV(:,:,i)=(h(i)/2)*M*((rx(i)*Dr([3,5,6],:)+sx(i)*Ds([3,5,6],:))*n(i,1)+(ry(i)*Dr([3,5,6],:)+sy(i)*Ds([3,5,6],:))*n(i,2));
end

UGV=zeros(size(M2D,1),size(M,2),size(rx,1));
for i=1:size(rx,1)
    UGV(:,:,i)=GUV(:,:,i)';
end

UV=zeros(size(M,1),size(M,2),size(h,1));
for i=1:size(h,1)
    UV(:,:,i) = h(i)/2*M;
end
