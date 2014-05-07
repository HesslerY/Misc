clear
tic
K=[];
V=[];
Kfinal=[];
sigma=0.08;

Nexp=100000;
N=100;
vvv=-30:.5:1;
a=0
for f=10.^(vvv/10);
	a=a+1;
	     if mod(round(a/length(vvv)*100),10)==0
         	disp(a/length(vvv)*100)
     end
    v=f;
    

    
    S=(complex(sigma*randn(Nexp,N),sigma*randn(Nexp,N)+(v)));

    %function param_rice(R,frequ)
    Ri=S';
    KK=[];
    DSR=[];
    ssigma=[];
    vv=[];
    v=[];
    p=1;
    Z=[];
    theta=[];
    thetas=[];
    thetat=[];

    for k=1:1:length(Ri(1,:))

        [mu,si]=normfit(real(Ri(:,k)));
        [mm,ss]=normfit(imag(Ri(:,k)));
        temp=sqrt(mu^2+mm^2);

        ssigma=[ssigma;(si+ss)/2];
        vv=[vv;sqrt(mu^2+mm^2)];
        theta(k)=atan2(mm,mu);
        %Z=[Z;complex(mu,mm)];

    end

K=[K;f/ssigma/sqrt(2)];
V=[V vv];

end

 for i=1:1:length(V(1,:))
     moy=mean(V(:,i));
     sig=std(V(:,i));
     min95=prctile(V(:,i),2.5);
     max95=prctile(V(:,i),97.5);
  
    Kfinal=[Kfinal;moy sig moy/sigma/sqrt(2) min95/sigma/sqrt(2) max95/sigma/sqrt(2)];
 end
%load ResultN100-10000.mat
%  figure,
% plot(20*log10(K),20*log10(Vfinal(:,5)),20*log10(K),20*log10(Vfinal(:,6)))
figure,
hold on
%plot(20*log10(K/2),20*log10(Vfinal(:,5)),20*log10(K/2),20*log10(Vfinal(:,6)),20*log10(K/2),20*log10(Vfinal(:,1)/sigma),20*log10(K/2),20*log10(K/2))
plot(20*log10(K),20*log10(Kfinal(:,3)),'b.')
plot(20*log10(K),20*log10(Kfinal(:,4)),'g.')
plot(20*log10(K),20*log10(Kfinal(:,5)),'r.')
plot(20*log10(K),20*log10(K),'k')
%plot(20*log10(Result(:,1)),20*log10(Result(:,2)),'bo')
xlabel('K')
ylabel('K estime')
grid on
% 
% figure,
% plot(10*log10(K),Kfinal(:,4))

save K100b.mat 'Kfinal'
     



    toc
