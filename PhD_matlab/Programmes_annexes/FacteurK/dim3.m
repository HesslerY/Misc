clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
N=30; %position de brasseur
Nfreq=50; %nombre de fréquences
Nexp=1000; %nombre d'expérience de la simulation de MC
a=0;
K=-40:1:30;
f0=1e9;
v2=2*sigma^2*10.^(K/10);
a=0
for f=1:1:length(v2);
    a=a+1;
    if mod(round(a/length(v2)*500000),10)==0
        disp(a/length(v2)*100)
    end
    
    %Génération du S21
    Sr=sigma*(randn(N,Nfreq,Nexp))+(sqrt(v2(f)))*cos(45);
    Si=sigma*(randn(N,Nfreq,Nexp))+(sqrt(v2(f)))*sin(45);

    %moyenne sur les positions de brasseur
    mr=mean(Sr,1).^2;
    mi=mean(Si,1).^2;


    vv=(mr+mi); %v^2

    ss=(var(Sr,0,1)+var(Si,0,1))*N/(N-1); %2*sigma^2
    %ss=2*0.08^2; %sigma déterministe

    %moyenne sur les fréquences
    k=mean(vv./ss,2);
    %k=mean(vv./ss,2)-1/N;

    %moyenne sur les expériences (lissage de Monte Carlo)
    kmoy=mean(k,3);

    min95=prctile(k,5);
    max95=prctile(k,95);

    Kfinal=[Kfinal;kmoy min95 max95];
end


figure(1)
hold on
plot(K,K,'LineWidth',1,'Color','k')
plot(K,10*log10(Kfinal(:,1)),'b')
plot(K,10*log10(Kfinal(:,2)),'g')
plot(K,10*log10(Kfinal(:,3)),'g')
legend('True values','Estimated values','Range of estimated values')%,'Intervalle de confiance analytique')
title(['K-factor estimation for ',num2str(N),' stirrer positions, ',num2str(Nfreq),' frequencies and ',num2str(Nexp), ' Monte Carlo experiments'])
xlabel('K')
axis equal
ylabel('K estime')
grid on
hold off
toc


