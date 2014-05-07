clear all

D1=[];
D2=[];
R=10;
Dm=1;
DM=4;
f=1e9;
lambda=3e8/f;
thetamin=pi/2-atan(DM/R);
thetamax=pi/2-atan(Dm/R);
for i=1:100
    H=[];
    tic
    disp(i)
    for j=1:300
        
        h=j/100;
        
        
        
        H=[H;h/lambda];
        
        %thetamax=pi;
        
        n=5;
        L=3;
        xa=L*rand(1,n);
        ya=L*rand(1,n);
        za=L*rand(1,n)+h/lambda;
        Ia=complex(rand(1,n),rand(1,n));
        
        
        
%         n = randint(1,1,[0,7]);
%         if mod(n,2)==0
%             n = n+1;
%         end
%         dx=lambda/4/lambda;
%         xa = dx:dx:dx*n;%D*rand(n,1);
%         ya = lambda*ones(1,n)/lambda;%D*rand(n,1);
%         za = lambda*ones(1,n)/lambda+h/lambda;%D*rand(n,1);
%         [Ia, dpha] = dolph3(1, dx, n, randint(1,1,[1,20]));
        S=2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   Diag
        nt = 180;
        np = 360;
        theta = linspace(0, pi, nt);
        phi = linspace(0, 2*pi, np);
        [thetaa, phia] = meshgrid(theta, phi);
        frant = fres(thetaa, phia, n, xa, ya, za, Ia);
        faant = fantres(frant, thetaa, phia, S);
        d1ant = Directivite(faant, thetaa, phia);
        %tracediag(faant,thetaa, phia);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        D1(i,j)=d1ant;
        %images
        
        xb=xa;
        yb=ya;
        zb=-za;
        Ib=1*Ia;
        
        S = 2;
        
        x=[xa xb];
        y=[ya yb];
        z=[za zb];
        I=[Ia Ib];
        
        
        N=length(x);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   Diag
        nt = round(180*(thetamax-thetamin)/pi);
        np = 360;
        theta = linspace(thetamin,thetamax, nt);
        phi = linspace(0, 2*pi, np);
        [thetaa, phia] = meshgrid(theta, phi);
        frant = fres(thetaa, phia, N, x, y, z, I);
        faant = fantres(frant, thetaa, phia, S);
        d2ant = Directivite(faant, thetaa, phia);
        %tracediag(faant,thetaa, phia);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        D2(i,j)=d2ant;
        
    end
    toc
end

%plot(H,D1,H,D2)
plot(H,10*log10(nanmean(D2)./nanmean(D1)))