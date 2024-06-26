function [GUGV,GUV,UGV,UV] = compute_matrix(rx,ry,sx,sy,J,h,n,Dr,Ds,M2D,M,en,sign)

GUGV=zeros(size(M2D,1),size(M2D,2),size(rx,1));
for i=1:size(rx,1)
    GUGV(:,:,i)=J(i)*(rx(i)^2+ry(i)^2)*Dr'*M2D*Dr+(rx(i)*sx(i)+ry(i)*sy(i))*(Dr'*M2D*Ds+Ds'*M2D*Dr)+(sx(i)^2+sy(i)^2)*Ds'*M2D*Ds;
end

GUV=zeros(size(M,1),size(M2D,2),size(rx,1));
for i=1:size(rx,1)
    if sign(i)==1
        GUV(:,:,i)=(h(i)/2)*M*((rx(i)*Dr(en(1,:,1),:)+sx(i)*Ds(en(1,:,1),:))*n(i,1)+(ry(i)*Dr(en(1,:,1),:)+sy(i)*Ds(en(1,:,1),:))*n(i,2));
    else
        GUV(:,:,i)=flip((h(i)/2)*M*((rx(i)*Dr(en(1,:,1),:)+sx(i)*Ds(en(1,:,1),:))*n(i,1)+(ry(i)*Dr(en(1,:,1),:)+sy(i)*Ds(en(1,:,1),:))*n(i,2)));
    end
end

UGV=zeros(size(M2D,1),size(M,2),size(rx,1));
for i=1:size(rx,1)
    UGV(:,:,i)=GUV(:,:,i)';
end

UV=zeros(size(M,1),size(M,2),size(h,1));
for i=1:size(h,1)
    UV(:,:,i) = h(i)/2*M;
end
