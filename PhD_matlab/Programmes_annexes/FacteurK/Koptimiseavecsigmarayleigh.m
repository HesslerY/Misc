clear
tic
K=[];
V2=[];
Kfinal=[];

Nexp=1000000;
N=50;
SIGMA=[];
K=-40:1:30
f0=1e9;
sigma=sqrt(2.13e20/2/(f0)^2.5)
B=sigma*sqrt(2/pi);
v2=2*sigma^2*10.^(K/10);

%2 possibiltés laquelle choisir ?:
%1ere
for i=1:1:Nexp
    SIGMA=[SIGMA;mean(raylrnd(B,N,1))*ones(1,N)];
end
% 2eme        
% SIGMA=raylrnd(B,Nexp,N);%*ones(1,N)];

a=0

for f=1:1:length(v2);
	a=a+1;
	     if mod(round(a/length(v2)*100),10)==0
         	disp(a/length(v2)*100)
     end

    Sr=SIGMA.*(randn(Nexp,N));
    Si=SIGMA.*(randn(Nexp,N))+(sqrt(v2(f)));

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

figure(1)
hold on
plot(K,K,'LineWidth',1,'Color','k')
plot(K,10*log10(Kfinal(:,3)),'b')
plot(K,10*log10(Kfinal(:,4)),'--g')
plot(K,10*log10(Kfinal(:,5)),'--r')

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
