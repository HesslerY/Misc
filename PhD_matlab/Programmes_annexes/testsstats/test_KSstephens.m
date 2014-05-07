%function [RR]=test_KS(S21_2f,frequencesf);
%S=[0.35975571306374,0.111496146323539,0.187472450205357,0.208104312809706,0.17404532437845,0.186083415018642,0.129096185121792,0.300367377805913,0.310325430401699,0.127274556550003,0.207031281300677,0.305972732741008,0.22771434781322,0.280227311133301,0.299472404448223,0.30309801682789,0.140473369041965,0.245589528280422,0.349665263143195,0.287059894072648,0.36668283557456,0.306383454881297,0.137592590036673,0.244121630831846,0.230163083408265,0.178462507168873,0.170621435376684,0.1489091434466,0.259862944070523,0.265494501882807,0.27398476079702,0.352845514015695,0.243139861036811,0.402463290171911,0.291668649724306,0.169438863467033,0.258194029018876,0.305001656587304,0.360960429191345,0.277207479576219,0.230922495762106,0.313626163010996,0.304325110493696,0.207605315021076,0.294332086264818,0.303767090852515,0.137246901185418,0.270781303845003,0.285166542757035,0.187811039805971;]

S=abs(S21_2);
%Test de KS sur le nombre de positions de brasseur (N) pour loi expo et
%Rayleigh
 %position de brasseur ou nombre d'echantillons
SN=[];

%POURRI=[]; %matrice qui r�cup�re les fr�quences ou le test �choue
D=[];%matrice des valeurs de d
M=0;
for f=1:1:length(S(:,1)) %Boucle sur les fr�quences
    R=abs(S(f,:)); % on prend le module des valeurs de la fr�quence consid�r�e
    
    
%     %--------------------------------------------------%
%     %Partie utile pour la sonde de champ (doublons....)
%     doublons=[];
%     for f=1:1:length(S(:,1))
%         R=abs(S(f,:));
%         R=sort(R);
%         for i=1:1:length(R)-1
%             if R(i)==R(i+1)
%                 doublons=[doublons,i];
%             end
%         end
%         if length(doublons)~=0
%             for j=1:1:length(doublons)
%                 R(doublons(j))=[];
%             end
%         end
%     end
%     %--------------------------------------------------%



    b=raylfit(R); %on r�cup�re les param�tres b
    %x= min(R):(max(R)-min(R))/(length(R)-1):max(R); %on g�n�re le tableau de valeurs
    y=sort(R);
    p = raylcdf(y,b); %cdf th�orique pour les valeurs g�n�r�es

    %cdf de la mesure
    
        SN=1/(length(R)):1/(length(R)):1;

    D=abs(p-SN);
    d=max(D); %calcul de l'ecart maximal
    N=length(R);
    if ((d-0.2/N)*(sqrt(N)+0.26+0.5/sqrt(N)))<1.094 %test, valeurs critiques de stephens
        M=M+1; %increment
    else
        %        POURRI=[POURRI;frequencesf(f)];
    end
end

RR=1-M/length(S(:,1)) %proportion de r�ussite