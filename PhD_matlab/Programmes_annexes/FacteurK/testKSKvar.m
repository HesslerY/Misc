clear
tic
K=[];
V2=[];
Nexp=1000;%nombre d'expérience de la simulation de MC

sigma=.08;
N=50; %position de brasseur

a=0;
K=-40:1:30;
f0=1e9;
v2=2*sigma^2*10.^(K/10);
a=0;
disp('K    RR')
for ff=1:1:length(v2);
    a=a+1;
    %Génération du S21
    Sr=sigma*(randn(N,Nexp))+(sqrt(v2(ff)))*cos(45);
    Si=sigma*(randn(N,Nexp))+(sqrt(v2(ff)))*sin(45);
    S=abs(complex(Sr,Si));

    SN=[];

    D=[];%matrice des valeurs de d
    M=0;
    for f=1:1:length(S(:,1)) %Boucle sur les fr�quences
        R=abs(S(f,:)); % on prend le module des valeurs de la fr�quence consid�r�e

        b=raylfit(R); %on r�cup�re les param�tres b
        y=sort(R);
        p = raylcdf(y,b); %cdf th�orique pour les valeurs g�n�r�es


        SN=1/(length(R)):1/(length(R)):1;

        D=abs(p-SN);
        d=max(D); %calcul de l'ecart maximal
        N=length(R);
        if ((d-0.2/N)*(sqrt(N)+0.26+0.5/sqrt(N)))<1.094 %test, valeurs critiques de stephens
            M=M+1; %increment
        end
    end

    RR(ff)=1-M/length(S(:,1));%proportion de r�ussite
disp([K(ff) RR(ff)])
end

figure(1)
plot(K,RR)
xlabel('K')
ylabel('Taux de rejet')
