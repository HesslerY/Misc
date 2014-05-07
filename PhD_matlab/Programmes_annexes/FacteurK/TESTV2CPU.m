clear
%GPUstart
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
Nexp=10000;
N=30;

a=0
K=-40:1:30

v2=2*sigma^2*10.^(K/10);
a=0
for f=1:1:length(v2);
	a=a+1;
	     if mod(round(a/length(v2)*100),10)==0
         	disp(a/length(v2)*100)
     end

    Sr=sigma*(randn(Nexp,N));
    Si=sigma*(randn(Nexp,N))+(sqrt(v2(f)));

    %function param_rice(R,frequ)
    Rr=Sr';
    Ri=Si';
    
    
    ssigma=[];
    vv=[];
    

        mu=sum(Rr)/N;
        %si=std(Rr);
        mm=sum(Ri)/N;
        %ss=std(Ri);

        %ssigma=(si+ss)/2;
        vv=(mu.^2+mm.^2);

     moy=mean(vv);
     sig=std(vv);
     min95=prctile(vv,2.5);
     max95=prctile(vv,97.5);
  
    Kfinal=[Kfinal;moy sig/sigma^2/2 moy/sigma^2/2 min95/sigma^2/2 max95/sigma^2/2];
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
