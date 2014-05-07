clear
disp('POS')
tic
oo=1
freq=[];
VC=0.757
MM=150; %nombre d'???chantillon
l=3.7; %longueur du segment
p=3.7;
h=2.9;
c=3e8;%
Lt=3e-6;
dmax=c*Lt
ordre=round(dmax/min([l;p;h]))


filename = sprintf('%delem_%dnscouloir111.mat',1,round(Lt/(1e-9)));

load(filename)



length(POS);
set(gca,'nextplot','replacechildren');
%RR=[0.8;0.90;0.92;0.95;0.98;0.99;0.995;0.999;1]
%for  uuu=1:1:length(RR)
disp('RI')
tic
Rx=0.998;
Ry=Rx;
Rz=Rx;
%ordre=100 %ordre de calcul... nombre de r?flexions envisag?es
N=20000*3
d=1; %dirac
for ee=1:MM
    X_1=l/2*(1+.3*randn);
    Y_1=p/2*(1+.3*randn);
    Z_1=h/2*(1+.3*randn);
    %position du r???cepteur
    %dmax=c*Lt
    %ordre=round(dmax/min([l;p;h]))
    
    Sx=zeros(1,N);
    Sy=zeros(1,N);
    Sz=zeros(1,N);
    for j=1:1:length(POS)
        
        DX=X_1-POS(j,1);%D(j,1);
        DY=Y_1-POS(j,2);%D(j,2);
        DZ=Z_1-POS(j,3);%D(j,3);
        
        dist=sqrt(DX^2+DY^2+DZ^2);
        zl=round((N-1)*dist/c/Lt);
        
        if zl<N
            alpha=POS(j,9);
            beta=POS(j,10);
            ca=cos(alpha);
            sa=sin(alpha);
            cb=cos(beta);
            sb=sin(beta);
            %               Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
            %                         -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
            %                         -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];
            %
            %
            distx=((-sb)^2+(1-(-sb)^2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ;
            disty=(-sb*cb*(1-ca))*DX+((cb)^2+(1-cb^2)*ca)*DY+(sb*sa)*DZ;
            distz=(-cb*sa)*DX+(-sb*sa)*DY+(ca)*DZ;
            
            %             POSreperesource=Ralpbeta*[X_1-POS(j,1);Y_1-POS(j,2);Z_1-POS(j,3)];
            %             distx=(POSreperesource(1));
            %             disty=(POSreperesource(2));
            %             distz=(POSreperesource(3));
            
            
            
            
            distxy=sqrt(distx^2+disty^2);
            %o=POS(j,4)+POS(j,5)+POS(j,6);
            R=Rx^POS(j,4)*Ry^POS(j,5)*Rz^POS(j,6);
            E=POS(j,8)*d/dist;%*(4*o^2+2)^(.5);
            costheta=distz/dist;
            sintheta=distxy/dist;
            cosphi=distx/distxy;
            sinphi=disty/distxy;
            Antth=-sintheta;
            
            %             Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
            %                         (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
            %                         (cb*sa) (sb*sa) ca];
            
            %V=Ralpbetainv*[(-sintheta*costheta*cosphi);(-sintheta*costheta*sinphi);(sintheta^2)]; %C'est ici !
            
            Vx=(((-sb)^2+(1-(-sb)^2)*ca)*(Antth*costheta*cosphi)+(-sb*cb*(1-ca))*(Antth*costheta*sinphi)+(-cb*sa)*(-sintheta*Antth));
            Vy=((-sb*cb*(1-ca))*(Antth*costheta*cosphi)+((cb)^2+(1-(cb)^2)*ca)*(Antth*costheta*sinphi)+(-sb*sa)*(-sintheta*Antth));
            Vz=((cb*sa)*(Antth*costheta*cosphi)+(sb*sa)*(Antth*costheta*sinphi)+ca*(-sintheta*Antth));
            
            Sx(zl+1)=Sx(zl+1)+R*E*Vx;
            Sy(zl+1)=Sy(zl+1)+R*E*Vy;
            Sz(zl+1)=Sz(zl+1)+R*E*Vz;
            
        end
    end    %toc
    %Signal=S;%conv(S,[0 0 0 .05 .1 .3 1 1 .3 .1 .05 0 0 0]);
    %---------FFT-----------------------------%
    %disp('FFT 1')
    Fs = N/Lt;                    % Sampling frequency
    T = 1/Fs;                     % Sample time
    L = N; %400000;                     % Length of signal
    tt = (0:L-1)*T;                % Time vector
    
    NFFT = 2^nextpow2(N); % Next power of 2 from length of y
    Yx = fft((Sx),NFFT)/N;
    Yy = fft((Sy),NFFT)/N;
    Yz = fft((Sz),NFFT)/N;
    f = Fs/2*linspace(0,1,NFFT/2);
    Ax=abs(Yx(1:NFFT/2));
    Ay=abs(Yy(1:NFFT/2));
    Az=abs(Yz(1:NFFT/2));
    freq(:,oo)=f;
    SSx(:,ee)=Ax;
    SSy(:,ee)=Ay;
    SSz(:,ee)=Az;
end
%Test
toc

disp('TESTS X')
tic

M=0;
for f=1:1:length(SSx(:,1)) %Boucle sur les fr?quences
    Rxx=(SSx(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    %b=raylfit(Rxx); %on r?cup?re les param?tres b
    parmhat = wblfit(Rxx);
    x= min(Rxx):(max(Rxx)-min(Rxx))/(length(Rxx)-1):max(Rxx); %on g?n?re le tableau de valeurs
    y=sort(Rxx)';
    ppx = wblcdf(y,parmhat(1),parmhat(2)); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppx)+log(1-flipud(ppx));
    A=[];
    for i=1:length(Rxx)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinalx(f,oo)=1;
    else
        RRfinalx(f,oo)=0;
    end
end
disp('TESTS Y')

M=0;
for f=1:1:length(SSy(:,1)) %Boucle sur les fr?quences
    Ryy=(SSy(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    parmhat = wblfit(Ryy); %on r?cup?re les param?tres b
    x= min(Ryy):(max(Ryy)-min(Ryy))/(length(Ryy)-1):max(Ryy); %on g?n?re le tableau de valeurs
    y=sort(Ryy)';
    ppy = wblcdf(y,parmhat(1),parmhat(2)); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppy)+log(1-flipud(ppy));
    A=[];
    for i=1:length(Ryy)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinaly(f,oo)=1;
    else
        RRfinaly(f,oo)=0;
    end
end
disp('TESTS Z')

RR1=[];
M=0;
for f=1:1:length(SSz(:,1)) %Boucle sur les fr?quences
    Rzz=(SSz(f,:)); % on prend le module des valeurs de la fr?quence consid?r?e
    
    
    parmhat = wblfit(Rzz); %on r?cup?re les param?tres b
    x= min(Rzz):(max(Rzz)-min(Rzz))/(length(Rzz)-1):max(Rzz); %on g?n?re le tableau de valeurs
    y=sort(Rzz)';
    ppz = wblcdf(y,parmhat(1),parmhat(2)); %cdf th?orique pour les valeurs g?n?r?es
    X=log(ppz)+log(1-flipud(ppz));
    A=[];
    for i=1:length(Rzz)
        A=[A;(2*i-1)*X(i)];
    end
    A2=-sum(A)/MM-MM;
    if A2*(1+0.2/sqrt(MM))>VC %test, valeurs critiques de stephens
        M=M+1; %increment
        RRfinalz(f,oo)=1;
    else
        RRfinalz(f,oo)=0;
    end
end


b=50;
Dx=[];
Dy=[];
Dz=[];
for i=1:1:length(RRfinalx)-b
    Dx(i)=sum(RRfinalx(i:i+b))/b;
    Dy(i)=sum(RRfinaly(i:i+b))/b;
    Dz(i)=sum(RRfinalz(i:i+b))/b;
    
end
toc
save('weibullcouloir.mat','freq','RRfinalx','RRfinaly','RRfinalz','SSx','SSy','SSz')

figure(1)
% hold on
plot(freq,cumsum(RRfinalx),freq,cumsum(RRfinaly),freq,cumsum(RRfinalz))%,freq,m2*freq+p2,'k',freq,m1*freq+p1,'k')
ylabel('rejets cumul???s')
xlabel('Fr???quence')
title([num2str(sum(RRfinalx+RRfinaly+RRfinalz)),' rejets pour ',num2str(3*length(freq)),' tests, ',num2str(MM),' ???chantillons,','R=',num2str(Rx),', ',num2str(l),'x',num2str(p),'x',num2str(h)])
grid on
xlim([0 max(freq)])
