function DG_plot(i,j)
M=2.^(1:i);
N=2.^(0:j);
%M=2:i+1;
%N=1:j+1;
tiledlayout(i,j+1);
for k=1:i
    for l=1:j+1
        nexttile
        [~,~,~,~,u,c4n2] = DG2(M(k),N(l));
        %[~,~,~,~,u,~,~,c4n2] = DG(0,1,0,1,M(k),M(k),N(l));

        plot3(c4n2(:,1),c4n2(:,2),u,'.')
    end
end

