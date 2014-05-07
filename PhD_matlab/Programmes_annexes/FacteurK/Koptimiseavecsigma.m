clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
Nexp=1000000;
N=30;

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

    Sr=sigma*(randn(Nexp,N));
    Si=sigma*(randn(Nexp,N))+(sqrt(v2(f)));

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
    min95simple=prctile(vv./(2*sigma^2),2.5);
    max95simple=prctile(vv./(2*sigma^2),97.5);
    min95=prctile(vv./(2*s.^2),2.5);
    max95=prctile(vv./(2*s.^2),97.5);
    
    Kfinal=[Kfinal;moy sig/sigmamoy^2/2 moy/sigmamoy^2/2 min95 max95 min95simple max95simple];
 end

figure(1)
hold on
plot(K,K,'LineWidth',1,'Color','k')
plot(K,10*log10(Kfinal(:,3)),'b')
plot(K,10*log10(Kfinal(:,4)),'--g')
plot(K,10*log10(Kfinal(:,5)),'--r')
% plot(K,10*log10(Kfinal(:,6)),'--g')
% plot(K,10*log10(Kfinal(:,7)),'--r')
plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(30)),'--m')
plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(100)),'--g')
plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(30)),'--m')
plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(100)),'--g')
legend('True values','Estimated values','Centile @ 97.5%','Centile @ 2.5%','Range of estimated values with N=30','Range of estimated values with N=100')
xlabel('K')
axis equal
ylabel('K estime')
grid on
toc