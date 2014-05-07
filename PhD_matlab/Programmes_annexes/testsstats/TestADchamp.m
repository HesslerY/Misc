
VC=1.341
MM=100; %nombre d'�chantillon/pos de brasseur
%SSx c'est le champ sur MM positions (colonnes) pour un nombre de fréquence (lignes)

disp('TESTS X')
tic

M=0;
for f=1:1:length(SSx(:,1)) %Boucle sur les fr?quences
    Rxx=abs(SSx(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    b=raylfit(Rxx); %on r?cup?re les param?tres b
    x= min(Rxx):(max(Rxx)-min(Rxx))/(length(Rxx)-1):max(Rxx); %on g?n?re le tableau de valeurs
    y=sort(Rxx)';
    ppx = raylcdf(y,b); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppx)+log(1-flipud(ppx));
    A=[];
    for i=1:length(Rxx)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.6/MM)>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinalx(f,oo)=1;
    else
        RRfinalx(f,oo)=0;
    end
end
disp('TESTS Y')

M=0;
for f=1:1:length(SSy(:,1)) %Boucle sur les fr?quences
    Ryy=abs(SSy(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    b=raylfit(Ryy); %on r?cup?re les param?tres b
    x= min(Ryy):(max(Ryy)-min(Ryy))/(length(Ryy)-1):max(Ryy); %on g?n?re le tableau de valeurs
    y=sort(Ryy)';
    ppy = raylcdf(y,b); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppy)+log(1-flipud(ppy));
    A=[];
    for i=1:length(Ryy)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.6/MM)>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinaly(f,oo)=1;
    else
        RRfinaly(f,oo)=0;
    end
end
disp('TESTS Z')

RR1=[];
M=0;
for f=1:1:length(SSz(:,1)) %Boucle sur les fr?quences
    Rzz=abs(SSz(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    b=raylfit(Rzz); %on r?cup?re les param?tres b
    x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
    y=sort(Rzz)';
    ppz = raylcdf(y,b); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppz)+log(1-flipud(ppz));
    A=[];
    for i=1:length(Rzz)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.6/MM)>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinalz(f,oo)=1;
    else
        RRfinalz(f,oo)=0;
    end
end


b=50;
Dx=[];
Dy=[];
Dz=[];
for i=1:1:length(RRfinalx)-b
    Dx(i)=sum(RRfinalx(i:i+b))/b;
    Dy(i)=sum(RRfinaly(i:i+b))/b;
    Dz(i)=sum(RRfinalz(i:i+b))/b;
    
end
toc
save('10011109983mu.mat','freq','RRfinalx','RRfinaly','RRfinalz')

figure(1)
% hold on
plot(freq,cumsum(RRfinalx),freq,cumsum(RRfinaly),freq,cumsum(RRfinalz))%,freq,m2*freq+p2,'k',freq,m1*freq+p1,'k')
ylabel('rejets cumul�s')
xlabel('Fr�quence')
title([num2str(sum(RRfinalx+RRfinaly+RRfinalz)),' rejets pour ',num2str(3*length(freq)),' tests, ',num2str(MM),' �chantillons,','R=',num2str(Rx),', ',num2str(l),'x',num2str(p),'x',num2str(h)])
grid on
xlim([0 max(freq)])
