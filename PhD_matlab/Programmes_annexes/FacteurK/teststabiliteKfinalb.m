clear
tic
K=[];
V=[];
S=[];
Kfinal=[];
sigma=0.08;

Nexp=10000;
N=100;
KdB=-40:.5:20;
a=0
for f=10.^(KdB/20)*sigma*sqrt(2);
	a=a+1;
	     if mod(round(a/length(KdB)*100),10)==0
         	disp(a/length(KdB)*100)
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
        
        ssigma=[ssigma;(si+ss)/2];
        vv=[vv;sqrt(mu^2+mm^2)];
        theta(k)=atan2(mm,mu);


    end


V=[V vv];
S=[S ssigma];
end

 for i=1:1:length(V(1,:))
     moy=mean(V(:,i));
     sig=std(V(:,i));

  
    Kfinal=[Kfinal;20*log10(moy/sqrt(2)/sigma)];
 end

figure,
hold on

plot(KdB,KdB,'k')
plot(KdB,Kfinal,'r.')

xlabel('K')
ylabel('K estime')
grid on


save K100b.mat 'Kfinal'
     



    toc
