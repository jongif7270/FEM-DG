function DGerror_plot(xl,xr,yl,yr,M,N)
K=2.^(1:M);
tiledlayout(M,1);
for j=1:M
    nexttile
    DGerror2(xl,xr,yl,yr,K(j),N);
end
