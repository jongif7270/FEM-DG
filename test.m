function [answer2,answer3,answer4,inddb2,inddb3,inddb4]=test()
M=2;N=2;
xl=-1;xr=1;yl=-1;yr=1;Mx=M;My=M;

[c4n,n4e,~,~] = mesh_fem_2d_triangle(xl,xr,yl,yr,Mx,My,N);
[ind4e,~,~,c4n2,~] = indexforDG(xl,xr,yl,yr,Mx,My,N);

ind4s=edge(n4e,N);
e4s = computeE4s(n4e);
s4e = computeS4e(n4e);

count=sum(e4s(:, 2) == 0);
inddb2=zeros(1,count*(N+1));
tmp=0;
for i=1:size(e4s,1)
    if e4s(i,2)==0
        tmp=tmp+1;
        inddb2((tmp-1)*(N+1)+(1:N+1))=(i-1)*(N+1)+(1:N+1);
    end
end

inddb3=[];
inddb4=[];
for i=1:size(e4s,1)
    if e4s(i,2)==0
        if mod(e4s(i,1),2)~=0
            inddb3=[inddb3 (i-1)*(N+1)+(1:N+1)];
        else
            inddb4=[inddb4 (i-1)*(N+1)+(1:N+1)];
        end
    end
end

answer2=zeros(size(e4s,1)*(N+1));

for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        if e4s(s4e(j,i),2)==0
            answer2((s4e(j,i)-1)*(N+1)+(1:N+1))=(s4e(j,i)-1)*(N+1)+(1:N+1);
        end
    end
end

answer2=nonzeros(answer2);

answer3=zeros(size(e4s,1)*(N+1));
answer4=zeros(size(e4s,1)*(N+1));
a=[0.8,0.6];
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        n=normal(c4n(n4e(j,mod(i,3)+1),:)-c4n(n4e(j,mod(i-1,3)+1),:));
        if e4s(s4e(j,i),2)==0
            if dot(a,n)>0
                answer4((s4e(j,i)-1)*(N+1)+(1:N+1))=(s4e(j,i)-1)*(N+1)+(1:N+1);
            else
                answer3((s4e(j,i)-1)*(N+1)+(1:N+1))=(s4e(j,i)-1)*(N+1)+(1:N+1);
            end
            %% ex 3.4용
            %v=dot(c4n2(ind4s(s4e(j,i),:,1),:)',repmat(n,1,N+1));
            %vp=0;
            %vn=0;
            %for k=1:N+1
            %    if v(k)>0
            %        vp=vp+1;
            %    else
            %        vn=vn+1;
            %    end
            %end
            %if vp>vn
            %    answer4=[answer4 (s4e(j,i)-1)*(N+1)+(1:N+1)];
            %else
            %    answer3=[answer3 (s4e(j,i)-1)*(N+1)+(1:N+1)];
            %end
        end
    end
end

answer2=nonzeros(answer2)';
answer3=nonzeros(answer3)';
answer4=nonzeros(answer4)';


a=@(x) [exp(x(:,1))*(x(:,2)*cos(x(:,2))+sin(x(:,2))),-exp(x(:,1))*x(:,2)*sin(x(:,2))];

ae=zeros((N+1)*(N+2)/2,2,size(n4e,1));
for j=1:size(n4e,1)
    for k=1:(N+1)*(N+2)/2
        ae(k,:,j)=a(c4n2(((N+1)*(N+2)/2)*(j-1)+k,:));
    end
end

al=zeros(N+1,2,size(e4s,1));
for j=1:size(n4e,1)
    for i=1:size(n4e,2)
        for k=1:N+1
            al(k,:,s4e(j,i))=a(c4n2(ind4s(s4e(j,i),k,1),:));
        end
    end
end
