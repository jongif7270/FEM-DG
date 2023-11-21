function DGerror_plot(xl,xr,yl,yr,N,t)
K=2.^(0:N);
tiledlayout(t,1);
for j=1:N+1
    nexttile
    DGerror(xl,xr,yl,yr,K(j),t);
end
