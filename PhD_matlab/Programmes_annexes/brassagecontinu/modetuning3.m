%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    MODE TUNING  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

load SRopldiff2.mat
SR=SR';
for uuu=1:1:500%length(SR(1,:))
    tic
    disp(uuu)
    SRi=((SR(:,uuu)./mean(SR(:,uuu))).^2);
    ppp=0;
    for ooo=0:0.1:10
        ppp=ppp+1;
        critere=ooo; %en dBm
        M=100;
        for N=1:M %positions de brasseurs
            Crit(N,ppp,uuu)=(critere);
            powerdB=10*log10(SRi)+3;
            essairand=0;
            posrand=0;
            nbressai=100;
            
            for u=1:nbressai
                posirand=sort(randi(length(SRi),N,1)); 
                P=powerdB(posirand);
                pmax(N,uuu,ppp)=max(P);
                U=find(P>critere);
                if length(U)>0
                    essairand=essairand+1;
                    posrand=posrand+posirand(U(1))/length(SRi)*360;
                end
            end
                       
            NN(N,uuu,ppp)=N;
            effrand(N,uuu,ppp)=essairand/nbressai*100;
            tempsrand(N,uuu,ppp)=posrand/nbressai; %en degrés
            
        end
    end
    toc
end
save('resmtv90expo_1M_100b.mat','NN','effrand','tempsrand','lamb','pmax','Crit')