clear
tic


Nexp=1000;
set(gcf,'nextplot','replacechildren');
for N=30;
K=[];
V2=[];
Kfinal=[];
K=-4:1:47
f0=1e9;
sigma=sqrt(2.13e20/2/(f0)^2.5)
B=sigma*sqrt(2/pi);
v2=2*sigma^2*10.^(K/10);


a=0

for f=1:1:length(v2);
	a=a+1;
	     if mod(round(a/length(v2)*100),10)==0
         	disp(a/length(v2)*100)
         end

     reverbr=sigma*randn(Nexp,N);
     reverbi=sigma*randn(Nexp,N);   
         
    Sr=reverbr;%;.*(randn(Nexp,N));
    Si=reverbi+sqrt(v2(f));

    Rr=Sr';
    Ri=Si';
    
    
    ssigma=[];
    vv=[];
    
mr=mean(Rr);
mi=mean(Ri);
sr=std(Rr);
si=std(Ri);

vv=(mr.^2+mi.^2);
s=(sr+si)/2;
    sigmamoy=mean(s);
    moy=mean(vv);
    sig=std(vv);
    min95=prctile(vv./(2*s.^2),2.5);
    max95=prctile(vv./(2*s.^2),97.5);
    
    Kfinal=[Kfinal;moy sig/sigmamoy^2/2 moy/sigmamoy^2/2 min95 max95];
 end

figure(2)


plot(K,K,K,10*log10(Kfinal(:,3)),K,10*log10(Kfinal(:,3)-1/N),K,10*log10(Kfinal(:,4)),K,10*log10(Kfinal(:,5)))%,K,10*log10(Kfinal(:,5)))
% plot(K,K,'LineWidth',1,'Color','k')
% plot(K,10*log10(Kfinal(:,3)),'b')
% plot(K,10*log10(Kfinal(:,4)),'--g')
% plot(K,10*log10(Kfinal(:,5)),'--r')

% plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(30)),'--m')
% plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(100)),'--g')
% plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(30)),'--m')
% plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(100)),'--g')
legend('K','Kestime','Kcorr')%,'Range of estimated values with N=30','Range of estimated values with N=100')
xlabel('K')
axis equal
title(['N= ',num2str(N)])
ylabel('K estime')
% xlim([-40 30])
% ylim([-40 50])
grid on
getframe;
          filename = sprintf('corrN_%d.png',N);
  filename2 = sprintf('corrN_%d.eps',N);
      % saveas(gcf,filename)%
%       
     %saveas(gca,filename2,'epsc')

toc
end