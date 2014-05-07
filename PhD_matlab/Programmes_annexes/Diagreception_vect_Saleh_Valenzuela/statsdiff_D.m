clear all
MM=50
paramomni=[];
parama=[];
DdB=[];
for j=1:5
 filename=sprintf('FFT_50-6muantenneorientationdx_%d_elem.mat',j)
 load(filename)
 freq=f;


param6mux=[];
param6muy=[];
param6muz=[];
param6muxa=[];
param6muya=[];
param6muza=[];
FFTz=abs(FFTz);
FFTza=abs(FFTza);

for i=1:length(FFTx)
%     A=wblfit(abs(FFTx(:,i)));
%     param6mux(i)=A(2);
%     A=wblfit(abs(FFTy(:,i)));
%     param6muy(i)=A(2);
    A=wblfit(FFTz(:,i));
    param6muz(i)=A(2);
%     A=wblfit(abs(FFTxa(:,i)));
%     param6muxa(i)=A(2);
%     A=wblfit(abs(FFTya(:,i)));
%     param6muya(i)=A(2);
    A=wblfit(FFTza(:,i));
    param6muza(i)=A(2);
end
paramomni=[paramomni;param6muz];
parama=[parama;param6muza];
DdB=[DdB;10*log10(D)];
end
p=100;
q=1;
plot(freq(1:q:32768),filter(ones(p,1)/p,1,paramomni(:,1:q:32768)'),freq(1:q:32768),filter(ones(p,1)/p,1,parama(:,1:q:32768)'))

% 
% %teststatistiques
% %Rayleigh
% 
% b=100;
% VC=1.341;
% disp('TESTS Rayleigh')
% 
% M=0;
% for f=1:1:length(FFTz(1,:)) %Boucle sur les fr?quences
%     Rzz=(FFTz(:,f))'; % on prend le module des valeurs de la fr?quence consid?r?e
%     
%     
%     parmhat = raylfit(Rzz); %on r?cup?re les param?tres b
%     x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
%     y=sort(Rzz)';
%     ppz = raylcdf(y,parmhat); %cdf th?orique pour les valeurs g?n?r?es
%     X=log(ppz)+log(1-flipud(ppz));
%     A=[];
%     for i=1:length(Rzz)
%         A=[A;(2*i-1)*X(i)];
%     end
%     A2=-sum(A)/MM-MM;
%     if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
%         M=M+1; %increment
%         RRfinalz6mur(f)=1;
%     else
%         RRfinalz6mur(f)=0;
%     end
% end
% 
% 
% Dz6r=[];
% for i=1:1:length(RRfinalz6mur)-b
%     
%     Dz6r(i)=sum(RRfinalz6mur(i:i+b))/b;
%     
% end
% 
% M=0;
% for f=1:1:length(FFTza(1,:)) %Boucle sur les fr?quences
%     Rzz=(FFTza(:,f))'; % on prend le module des valeurs de la fr?quence consid?r?e
%     
%     
%     parmhat = raylfit(Rzz); %on r?cup?re les param?tres b
%     x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
%     y=sort(Rzz)';
%     ppz = raylcdf(y,parmhat); %cdf th?orique pour les valeurs g?n?r?es
%     X=log(ppz)+log(1-flipud(ppz));
%     A=[];
%     for i=1:length(Rzz)
%         A=[A;(2*i-1)*X(i)];
%     end
%     A2=-sum(A)/MM-MM;
%     if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
%         M=M+1; %increment
%         RRfinalz6mura(f)=1;
%     else
%         RRfinalz6mura(f)=0;
%     end
% end
% 
% 
% Dz6ra=[];
% for i=1:1:length(RRfinalz6mura)-b
%     
%     Dz6ra(i)=sum(RRfinalz6mura(i:i+b))/b;
%     
% end
% 
% 
% VC=0.757;
% disp('TESTS Weibull')
% 
% M=0;
% for f=1:1:length(FFTz(1,:)) %Boucle sur les fr?quences
%     Rzz=(FFTz(:,f))'; % on prend le module des valeurs de la fr?quence consid?r?e
%     parmhat = wblfit(Rzz); %on r?cup?re les param?tres b
%     x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
%     y=sort(Rzz)';
%     ppz = wblcdf(y,parmhat(1),parmhat(2)); %cdf th?orique pour les valeurs g?n?r?es
%     X=log(ppz)+log(1-flipud(ppz));
%     A=[];
%     for i=1:length(Rzz)
%         A=[A;(2*i-1)*X(i)];
%     end
%     A2=-sum(A)/MM-MM;
%     if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
%         M=M+1; %increment
%         RRfinalz6mu(f)=1;
%     else
%         RRfinalz6mu(f)=0;
%     end
% end
% 
% 
% Dz6=[];
% for i=1:1:length(RRfinalz6mu)-b
%     
%     Dz6(i)=sum(RRfinalz6mu(i:i+b))/b;
%     
% end
% 
% 
% M=0;
% for f=1:1:length(FFTza(1,:)) %Boucle sur les fr?quences
%     Rzz=(FFTza(:,f))'; % on prend le module des valeurs de la fr?quence consid?r?e
%     parmhat = wblfit(Rzz); %on r?cup?re les param?tres b
%     x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
%     y=sort(Rzz)';
%     ppz = wblcdf(y,parmhat(1),parmhat(2)); %cdf th?orique pour les valeurs g?n?r?es
%     X=log(ppz)+log(1-flipud(ppz));
%     A=[];
%     for i=1:length(Rzz)
%         A=[A;(2*i-1)*X(i)];
%     end
%     A2=-sum(A)/MM-MM;
%     if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
%         M=M+1; %increment
%         RRfinalz6mua(f)=1;
%     else
%         RRfinalz6mua(f)=0;
%     end
% end
% 
% 
% Dz6a=[];
% for i=1:1:length(RRfinalz6mua)-b
%     
%     Dz6a(i)=sum(RRfinalz6mua(i:i+b))/b;
%     
% end
% q=1
% 
% figure(1)
% plot(freq(1:4*q:8092)/1e6,Dz6r(1:4*q:8092),freq(1:4*q:8092)/1e6,Dz6ra(1:4*q:8092))
% legend('Omni','Antenne')
% xlabel('Frequency in MHz')
% ylabel('Reject rate in\%')
% grid on
% 
% figure(2)
% plot(freq(1:4*q:8092)/1e6,Dz6(1:4*q:8092),freq(1:4*q:8092)/1e6,Dz6a(1:4*q:8092))
% legend('Omni','Antenne')
% xlabel('Frequency in MHz')
% ylabel('Reject rate in\%')
% grid on
% 
% p=100
% 
% figure(3)
% plot(freq,filter(ones(p,1)/p,1,param6muz),freq,filter(ones(p,1)/p,1,param6muza))
% legend('Omni','Antenne')
% xlabel('Frequency in MHz')
% ylabel('$\beta$')
% grid on
% 
