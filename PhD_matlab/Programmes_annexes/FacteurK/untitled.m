clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
N=30; %position de brasseur
Nfreq=30; %nombre de fréquences
Nexp=10000; %nombre d'expérience de la simulation de MC
a=0;
K=-40:1:0;
f0=1e9;
v2=2*sigma^2*10.^(K/10);
a=0
%B=5:55
%A=[0.9:.01:1.1];
A=[28/29]
%for i=1:1:length(A)
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

        ss1=(var(Sr,1,1)+var(Si,1,1))/(N/(N-1));
        ss2=(var(Sr,1,1)+var(Si,1,1))/((N-1)/(N-2));
        %2*sigma^2
        %ss=2*0.08^2; %sigma déterministe

        %moyenne sur les fréquences
        k1=mean(vv./ss1,2)-1/N;
        k2=mean(vv./ss2,2)-1/N;
        %moyenne sur les expériences (lissage de Monte Carlo)
        kmoy1=mean(k1,3);
        kmoy2=mean(k2,3);


        Kfinal1(f)=kmoy1;
        Kfinal2(f)=kmoy2;
    end
%end

%save KfacteurAb.mat 'Kfinal'
%save Ab.mat 'A'
%save Kb.mat 'K'

%BBB=[A(14:2:25)';A(26)';A(28:2:40)';0]
figure(1)
hold on
plot(K,K,'k')

plot(K,10*log10(Kfinal1),'--b+')
plot(K,10*log10(Kfinal2),'--g+')
%plot(K,10*log10(Kfinal(1:1:21,:)),'--r')

legend('K','Nous','Nous2')%,'Autres entre 0.9 et 1.1')
%legend('K','20/19','22/21','24/23','26/25','28/27','30/29','32/31','34/33','36/35','38/37','40/39')
%title(['K-factor estimation for ',num2str(N),' stirrer positions, ',num2str(Nfreq),' frequencies and ',num2str(Nexp), ' Monte Carlo experiments'])
xlabel('K')
axis equal
ylabel('K estime')
grid on
hold off
toc
xlim([-40 0])
