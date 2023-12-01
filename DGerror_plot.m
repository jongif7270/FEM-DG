function DGerror_plot(xl,xr,yl,yr,t,N)
K=2.^(0:N);
tiledlayout(N+1,1);
for j=1:N+1
    nexttile
    DGerror(xl,xr,yl,yr,t,K(j));
end
