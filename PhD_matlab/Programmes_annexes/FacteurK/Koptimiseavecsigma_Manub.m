clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
Nexp=100000;
N=100;

a=0
K=-40:1:30
f0=1e9;
v2=2*sigma^2*10.^(K/10);
a=0
for f=1:1:length(v2);
	a=a+1;
	     if mod(round(a/length(v2)*100),10)==0
         	disp(a/length(v2)*100)
     end

    Rr=sigma*(randn(N,Nexp));
    Ri=sigma*(randn(N,Nexp))+(sqrt(v2(f)));

    
    
    
   
    vv=[];
    
mr=mean(Rr).^2;
mi=mean(Ri).^2;


vv=(mr+mi);
%s=(sr+si)/2;
ss=(var(Rr)+var(Ri));
%sigmamoy=mean(s);
 %   moy=mean(vv);
  %  sig=std(vv);
kmoy=mean(vv./ss);
  sig=std((vv./ss));
  
  
    min95=prctile(vv./ss,2.5);
    max95=prctile(vv./ss,97.5);
    
    Kfinal=[Kfinal;kmoy sig min95 max95];
end

%******************
% Kvrai_dB=-40:1:30;
% Kvrai_lin=10.^(Kvrai_dB/10);
% q025_lin=Kvrai_lin+1/N-1.96*sqrt(1/(N*N)+2*Kvrai_lin/N)/sqrt(N);
% q025_dB=10*log10(q025_lin);
% q975_lin=Kvrai_lin+1/N+1.96*sqrt(1/(N*N)+2*Kvrai_lin/N)/sqrt(N);
% q975_dB=10*log10(q975_lin);
% y=[Kvrai_dB;q025_dB;q975_dB];
%******************

figure(1)
hold on
plot(K,K,'LineWidth',1,'Color','k')
plot(K,10*log10(Kfinal(:,1)),'b')
 plot(K,10*log10(Kfinal(:,3)),'--g')
 plot(K,10*log10(Kfinal(:,4)),'-.r')
%plot(K,10*log10(Kfinal(:,1)-2*Kfinal(:,2)/sqrt(N)),'--m')
% plot(Kvrai_dB,y(2,:),'b.');
% plot(Kvrai_dB,y(3,:),'b.');
%plot(K,10*log10(Kfinal(:,1)-2*Kfinal(:,2)/sqrt(100)),'--g')
%plot(K,10*log10(Kfinal(:,1)+2*Kfinal(:,2)/sqrt(N)),'--m')
%plot(K,10*log10(Kfinal(:,1)+2*Kfinal(:,2)/sqrt(100)),'--g')
%plot(K,10*log10(10.^(K/10)+1/N),'b.')

legend('True values','Estimated values','Centile @ 97.5%','Centile @ 2.5%')
xlabel('K')
axis equal
ylabel('Estimated K')
grid on
hold off
toc
xlim([-40 30])
figure(2)
var_th=1/N^2+2*(10.^(K/10))/N;
plot(K,20*log10(Kfinal(:,2)),K,10*log10(var_th));

