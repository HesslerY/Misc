%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    MODE TUNING  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

load SR360_3000MHZ_360.mat
theta=1:360;

%Amplitude_crete=sin(theta/16.7).^4;
thetai = .001:.001:360;
Ampi=interp1(theta,Amplitude_crete,thetai,'spline'); 
%SR=(Ampi/mean(Ampi)).^2;
SR=(Ampi/mean(Ampi)).^2;


critere=4 %en dBm
ttt=NaN;
M=60;
for N=1:M %positions de brasseurs
    disp(N)
    TR=(ttt*(M-N)-sum((ttt/N)*(1:N)));
    tic
    
    
    disp(TR)
    
    
    
    powerdB=10*log10(SR);
    
    essairand=0;
    posrand=0;
    essaiequi=0;
    posequi=0;
    nbressai=100;
    
    for u=1:nbressai
        %         if mod(u,10000)==0
        %             disp(u)
        %         end
        
        %%% 1ère stratégie: 50 positions aléatoires
                posirand=sort(randi(length(SR),N,1)); %pas de positions differentes
%         posirand=[randi(length(SR))];
%         
%         while length(posirand)<N
%             a=randi(length(SR));
%             if min(abs(posirand-a))<1000
%                 posirand=[posirand;a];
%             end
%         end
        posirand=sort(posirand);
        for i=1:N;
            P(i)=powerdB(posirand(i));
        end
        
        [vmax,posmax]=max(P);
        if vmax>critere
            %disp('essai OK')
            essairand=essairand+1;
            %disp(posmax)
            posrand=posrand+posirand(posmax)/length(SR);
            %else
            % disp('echec')
        end
        
        %%% 2ème stratégie 50 positions équi-réparties
        
        START=randi(length(SR));
        for i=1:N;
            Q(i)=powerdB(mod(START+i*round((length(SR)/N)),length(SR))+1);
        end
        
        [vmax,posmax]=max(Q);
        if vmax>critere
            %     disp('essai OK')
            %     disp(posmax)
            essaiequi=essaiequi+1;
            
            posequi=posequi+posmax*round((length(SR)/N))/length(SR);
            %else
            %disp('echec')
        end
    end
    
    
    
    % disp('défaut rand')
    % essairand/nbressai*100
    % disp('temps rand en ° parcourus')
    % posrand
    %
    % disp('défaut equi')
    % essaiequi/nbressai*100
    % disp('temps rand en ° parcourus')
    % posequi
    
    NN(N)=N;
    effrand(N)=essairand/nbressai*100;
    effequi(N)=essaiequi/nbressai*100;
    tempsrand(N)=posrand;
    tempsequi(N)=posequi;
 ttt=toc;
end

figure(1)
plot(NN,effrand,'--k',NN,effequi,'k')
xlabel('$N$')
ylabel('Failure \%')
legend('Random positions','Uniform positions')
title('3GHz, $>$6dB')
grid on

figure(2)
plot(NN,tempsrand,'--k',NN,tempsequi,'k')
xlabel('$N$')
ylabel('Degrees before failure')
legend('Random positions','Uniform positions')
title('3GHz, $>$6dB')
grid on
