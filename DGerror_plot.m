function DGerror_plot(M,N)
K=2.^(0:N);
tiledlayout(N+1,1);
for j=1:N+1
    nexttile
    DGerror(M,K(j));
end
