clear all
load dire2.mat
load rand2.mat
plot(H,(nanmean(D2)./nanmean(D1)),H,(nanmean(D2r)./nanmean(D1r)))
legend('random but directive radiation pattern','random radiation pattern')

