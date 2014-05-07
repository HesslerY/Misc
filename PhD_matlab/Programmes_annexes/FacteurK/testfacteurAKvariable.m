clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
N=30; %position de brasseur
Nfreq=30; %nombre de fréquences
Nexp=100000; %nombre d'expérience de la simulation de MC
a=0;
K=-40:.5:0;
f0=1e9;
v2=2*sigma^2*10.^(K/10);
a=0
B=5:55
A=(B./(B-1));
for i=1:1:length(A)
for f=1:1:length(v2);
    a=a+1;
    if mod(round(a/length(A)/length(K)*500000),10)==0
        disp(a/length(A)/length(K)*100)
    end
    
    %Génération du S21
    Sr=sigma*(randn(N,Nfreq,Nexp))+(sqrt(v2(f)))*cos(45);
    Si=sigma*(randn(N,Nfreq,Nexp))+(sqrt(v2(f)))*sin(45);

    %moyenne sur les positions de brasseur
    mr=mean(Sr,1).^2;
    mi=mean(Si,1).^2;


    vv=(mr+mi); %v^2

    ss=(var(Sr,0,1)+var(Si,0,1))*A(i);

    %2*sigma^2
    %ss=2*0.08^2; %sigma déterministe

    %moyenne sur les fréquences
    k=mean(vv./ss,2)-1/N;

    %moyenne sur les expériences (lissage de Monte Carlo)
    kmoy=mean(k,3);



    Kfinal(i,f)=kmoy;
end
end

save KfacteurA.mat 'Kfinal'
save A.mat 'A'
save K.mat 'K'

BBB=[A(14:2:25)';A(26)';A(28:2:40)';0]
figure(1)
hold on
plot(K,K,'k')
plot(K,10*log10(Kfinal(16:2:25,:)),'--r')
plot(K,10*log10(Kfinal(26,:)),'--b+')
plot(K,10*log10(Kfinal(28:2:36,:)),'--r')



legend('K','20/19','22/21','24/23','26/25','28/27','30/29','32/31','34/33','36/35','38/37','40/39')
title(['K-factor estimation for ',num2str(N),' stirrer positions, ',num2str(Nfreq),' frequencies and ',num2str(Nexp), ' Monte Carlo experiments'])
xlabel('K')
axis equal
ylabel('K estime')
grid on
hold off
toc
xlim([-40 -20])

figure(2)
hold on
plot(K,K,'k')
plot(K,10*log10(Kfinal(1:3:51,:)),'r')
plot(K,10*log10(Kfinal(26,:)),'b')




%legend(num2str(A(1:3:51)'))
title(['K-factor estimation for ',num2str(N),' stirrer positions, ',num2str(Nfreq),' frequencies and ',num2str(Nexp), ' Monte Carlo experiments'])
xlabel('K')
axis equal
ylabel('K estime')
grid on
hold off
toc





