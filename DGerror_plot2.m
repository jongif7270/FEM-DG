function DGerror_plot2(M,t)
K=2.^(1:M);
tiledlayout(M,1);
for j=1:M
    nexttile
    DGerror(K(j),t);
end
