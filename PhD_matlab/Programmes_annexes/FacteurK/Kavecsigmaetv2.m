clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=0.08;
Nexp=1000;
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

    S=(complex(sigma*randn(Nexp,N),sigma*randn(Nexp,N)+(sqrt(v2(f)))));

    %function param_rice(R,frequ)
    Ri=S';
    
    ssigma=[];
    vv=[];
    

    for k=1:1:length(Ri(1,:))

        [mu,si]=normfit(real(Ri(:,k)));
        [mm,ss]=normfit(imag(Ri(:,k)));
       

        ssigma=[ssigma;(si+ss)/2];
        vv=[vv;(mu^2+mm^2)];
       % theta(k)=atan2(mm,mu);
        %Z=[Z;complex(mu,mm)];

    end

%K=[K;f/sigma/sqrt(2)];
V2=[V2 vv];

end

 for i=1:1:length(V2(1,:))
     moy=mean(V2(:,i));
     sig=std(V2(:,i));
     min95=prctile(V2(:,i),2.5);
     max95=prctile(V2(:,i),97.5);
  
    Kfinal=[Kfinal;moy sig/sigma^2/2 moy/sigma^2/2 min95/sigma^2/2 max95/sigma^2/2];
 end
%load ResultN100-10000.mat
%  figure,
% plot(20*log10(K),20*log10(Vfinal(:,5)),20*log10(K),20*log10(Vfinal(:,6)))
figure(1)
hold on
%plot(20*log10(K/2),20*log10(Vfinal(:,5)),20*log10(K/2),20*log10(Vfinal(:,6)),20*log10(K/2),20*log10(Vfinal(:,1)/sigma),20*log10(K/2),20*log10(K/2))
plot(K,K,'LineWidth',1,'Color','k')
plot(K,10*log10(Kfinal(:,3)),'b')
plot(K,10*log10(Kfinal(:,4)),'--g')
plot(K,10*log10(Kfinal(:,5)),'--r')
plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(30)),'--m')
plot(K,10*log10(Kfinal(:,3)-2*Kfinal(:,2)/sqrt(100)),'--g')
plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(30)),'--m')
plot(K,10*log10(Kfinal(:,3)+2*Kfinal(:,2)/sqrt(100)),'--g')
legend('True values','Estimated values','Centile @ 97.5%','Centile @ 2.5%','Range of estimated values with N=30','Range of estimated values with N=100')
%plot(20*log10(Result(:,1)),20*log10(Result(:,2)),'bo')
xlabel('K')
axis equal
ylabel('K estime')
grid on
% 
% figure,
% plot(10*log10(K),Kfinal(:,4))

%save TESTDUDSR1000.txt 'Vfinal'
     



    toc
