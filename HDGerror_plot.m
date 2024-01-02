function HDGerror_plot(M,N)
K=2.^(0:N);
tiledlayout(N+1,1);
for j=1:N+1
    nexttile
    HDGerror(M,K(j));
end
