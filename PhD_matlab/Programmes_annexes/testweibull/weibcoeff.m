parmhatx = [];
parmhaty = [];
parmhatz = [];
%load weibull1503mu.mat
for f=1:1:8000 %Boucle sur les fr?quences
    Rxx=(SSx(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    Ryy=(SSy(f,:));
    Rzz=(SSz(f,:));
    
    %b=raylfit(Rxx); %on r?cup?re les param?tres b
    parmhatx = [parmhatx;wblfit(Rxx)];
    parmhaty = [parmhaty;wblfit(Ryy)];
    parmhatz = [parmhatz;wblfit(Rzz)];
end
betachr=[200 1.46;300 1.56;400 1.53;500 1.68;600 1.69;700 1.8;800 1.8; 900 1.86;1000 1.88;1100 1.93]
p=40
bx=filter(ones(p,1)/p,1,parmhatx(:,2));
by=filter(ones(p,1)/p,1,parmhaty(:,2));
bz=filter(ones(p,1)/p,1,parmhatz(:,2));
q=40
figure(1)
plot(betachr(:,1),betachr(:,2),freq(1:q:8000)/1e6,bx(1:40:8000),freq(1:q:8000)/1e6,by(1:q:8000),freq(1:q:8000)/1e6,bz(1:q:8000))
grid on
ylabel('$\beta$')
xlabel('Frequency in MHz')
legend('Measurements','$E_x$','$E_y$','$E_z$')
xlim([0 1500])
ylim([0.8 3.5])
