clear
tic
K=[];
V2=[];
Kfinal=[];
sigma=.08;
Nexp=100000;

    K=-35%-40:1:30
    f0=1e9;
    v2=2*sigma^2*10.^(K/10);
N=10:300
for i=1:1:length(N)

disp(i/length(N)*100)
    for f=1:1:length(v2);

        Rr=sigma*(randn(N(i),Nexp));
        Ri=sigma*(randn(N(i),Nexp))+(sqrt(v2(f)));

        mr=mean(Rr).^2;
        mi=mean(Ri).^2;

        vv=(mr+mi);
        ss=(var(Rr)+var(Ri));
        
        kmoy=mean(vv./ss);
        sig=std((vv./ss));
        
        kbon=mean(vv)/mean(ss);
   
        min95=prctile(vv./ss,2.5);
        max95=prctile(vv./ss,97.5);

        Kfinal(i,f)=kbon/kmoy;
    end


end
theo=(N-1)./N;
figure(1)
hold on
plot(N(1:1:291),Kfinal(1:1:291),'.')
plot(N,theo,'r')
legend('Correction factor','(N-1)/N')
xlabel('N')
grid on

