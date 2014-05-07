%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    MODE TUNING  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

load SR360_3000MHZ_360.mat
theta=1:360;

%Amplitude_crete=sin(theta/16.7).^4;

SR=(Amplitude_crete/mean(Amplitude_crete)).^2;

critere=6 %en dBm

for N=1:100 %positions de brasseurs
    disp(N)
    powerdB=10*log10(SR);
    
    essairand=0;
    posrand=0;
    essaiequi=0;
    posequi=0;
    nbressai=100000;
    
    for u=1:nbressai
        %         if mod(u,10000)==0
        %             disp(u)
        %         end
        
        %%% 1ère stratégie: 50 positions aléatoires
        posirand=sort(randi(360,N,1));
        for i=1:N;
            P(i)=powerdB(posirand(i));
        end
        
        [vmax,posmax]=max(P);
        if vmax>critere
            %disp('essai OK')
            essairand=essairand+1;
            %disp(posmax)
            posrand=posrand+posirand(posmax)/360;
            %else
            % disp('echec')
        end
        
        %%% 2ème stratégie 50 positions équi-réparties
        
        START=randi(360);
        for i=1:N;
            Q(i)=powerdB(mod(START+i*round((360/N)),360)+1);
        end
        
        [vmax,posmax]=max(Q);
        if vmax>critere
            %     disp('essai OK')
            %     disp(posmax)
            essaiequi=essaiequi+1;
            
            posequi=posequi+posmax*round((360/N))/360;
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

