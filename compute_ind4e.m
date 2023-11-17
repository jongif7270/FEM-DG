function [ind4e] = compute_ind4e(c4n,n4e,k)
allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
[~,ind] = unique(sort(allSides,2),'rows','first');
number_of_sides=length(ind);
number_of_elements=length(n4e(:,1));
number_of_nodes=length(c4n(:,1));
ind4s=zeros(number_of_sides, k-1);
for j=1:k-1
    ind4s(:,j)=j:k-1:number_of_sides*(k-1);
end
ind4s=ind4s+number_of_nodes;
s4e = computeS4e(n4e);
tmp1=zeros(number_of_elements,(k-1));
tmp2=zeros(number_of_elements,(k-1));
tmp3=zeros(number_of_elements,(k-1));
for j=1:number_of_elements
    numbertocheck1=s4e(j,1);
    A=s4e(1:j,:);
    countofnumber1=sum(A(:)==numbertocheck1);
    if countofnumber1 == 2
        tmp1(j,:)=flip(ind4s(s4e(j,1),:));
    else 
        tmp1(j,:)=ind4s(s4e(j,1),:);
    end
    numbertocheck2=s4e(j,2);
    countofnumber2=sum(A(:)==numbertocheck2);
    if countofnumber2 == 2
        tmp2(j,:)=flip(ind4s(s4e(j,2),:));
    else 
        tmp2(j,:)=ind4s(s4e(j,2),:);
    end
    numbertocheck3=s4e(j,3);
    countofnumber3=sum(A(:)==numbertocheck3);
    if countofnumber3 == 2
        tmp3(j,:)=flip(ind4s(s4e(j,3),:));
    else 
        tmp3(j,:)=ind4s(s4e(j,3),:);
    end
end    
tmp=[tmp1 tmp2 tmp3];
ind4int=zeros(number_of_elements,(k-2)*(k-1)/2);
for j=1:(k-2)*(k-1)/2
    ind4int(:,j)=j:(k-2)*(k-1)/2:number_of_elements*((k-2)*(k-1)/2);
end
ind4int=ind4int+number_of_nodes+number_of_sides*(k-1);
ind4e=[n4e tmp ind4int];