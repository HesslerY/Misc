tic
K=[];
V=[];
Vfinal=[];
sigma=0.08;
for f=0:.01:4*sigma    
	a=(100*f/4);
	     if mod(a,5)==0
         	disp(a)
     end
    v=f*sigma;
    Nexp=100;
    N=50;

    
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
    %theta=angle(Z);
% 
%     figure,
%     plot(theta)
   
%     KK=(vv./(sqrt(2)*ssigma))
%     plot(20*log10(KK))
K=[K;f];
V=[V vv];



end

 for i=1:1:length(V(1,:))
     moy=mean(V(:,i));
     sig=std(V(:,i));
     min95=prctile(V(:,i),2.5);
     max95=prctile(V(:,i),97.5);
    % Vfinal=[Vfinal;moy 20*log10(moy/sig) sig sig/moy 20*log10((moy-2*sig)/sigma) 20*log10((moy+2*sig)/sigma)]; 
    Vfinal=[Vfinal;moy/2 (moy/sig/2) sig sig/moy ((min95)/sigma)/2 ((max95)/sigma/2)];
 end
load resstatN50-100b.mat
%  figure,
% plot(20*log10(K),20*log10(Vfinal(:,5)),20*log10(K),20*log10(Vfinal(:,6)))
figure,
hold on
%plot(20*log10(K/2),20*log10(Vfinal(:,5)),20*log10(K/2),20*log10(Vfinal(:,6)),20*log10(K/2),20*log10(Vfinal(:,1)/sigma),20*log10(K/2),20*log10(K/2))
plot(10*log10(K/2),10*log10(Vfinal(:,5)),'g.')
plot(10*log10(K/2),10*log10(Vfinal(:,6)),'r.')
plot(10*log10(K/2),10*log10(Vfinal(:,1)/sigma),'b.')
plot(10*log10(K/2),10*log10(K/2),'k')
plot(10*log10(Result(:,1)./0.16),10*log10(Result(:,2)./0.16),10*log10(Result(:,1)./0.16),10*log10(Result(:,1)./0.16))
xlabel('K')
ylabel('K estime')
grid on

figure,
plot(10*log10(K),Vfinal(:,4))

%save TESTDUDSR1000.txt 'Vfinal'
     



    toc
