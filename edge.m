function ind4s=edge(n4e,N)
e4s = computeE4s(n4e);
ind4s=zeros(size(e4s,1),N+1,2);
tmp=[];
for j=2:size(e4s,1)
    if e4s(j,1)<e4s(j-1,1)
        tmp=[tmp j-1];
    end
end

a=zeros(1,N+1);
a(1)=N+1;
for j=1:N
    a(j+1)=a(j)+N+1-j;
end

b=zeros(1,N+1);
b(1)=1;
for j=1:N
    b(j+1)=b(j)+N+2-j;
end

for j=1:tmp(1)
    ind4s(j,:,1)=ind4s(j,:,1)+(a+(e4s(j,1)-1)*(N+1)*(N+2)/2);
    if e4s(j,2)==0
        ind4s(j,:,2)=ind4s(j,:,1);
    else
        ind4s(j,:,2)=ind4s(j,:,2)+fliplr(a+(e4s(j,2)-1)*(N+1)*(N+2)/2);
    end
end

for j=tmp(1)+1:tmp(2)
    ind4s(j,:,1)=ind4s(j,:,1)+(b+(e4s(j,1)-1)*(N+1)*(N+2)/2);
    if e4s(j,2)==0
        ind4s(j,:,2)=ind4s(j,:,1);
    else
        ind4s(j,:,2)=ind4s(j,:,2)+fliplr(b+(e4s(j,2)-1)*(N+1)*(N+2)/2);
    end
end

for j=tmp(2)+1:size(e4s,1)
    ind4s(j,:,1)=ind4s(j,:,1)+((1:N+1)+(e4s(j,1)-1)*(N+1)*(N+2)/2);
    if e4s(j,2)==0
        ind4s(j,:,2)=ind4s(j,:,1);
    else
        ind4s(j,:,2)=ind4s(j,:,2)+fliplr(((1:N+1)+(e4s(j,1)-1)*(N+1)*(N+2)/2));
    end
end