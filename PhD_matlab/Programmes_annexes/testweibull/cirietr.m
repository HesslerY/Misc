%clear
disp('POS')
tic
oo=1
freq=[];
VC=0.757
MM=50; %nombre d'???chantillon

c=3e8;%
%Lt=2e-6;
%
%ordre=201;%round(dmax/min
Lt=min([l;p;h])*ordre/c;
dmax=c*Lt
%filename = sprintf('%delemIETR_%dns.mat',1,round(Lt/(1e-9)));

%load(filename)



length(POS);
set(gca,'nextplot','replacechildren');
%RR=[0.8;0.90;0.92;0.95;0.98;0.99;0.995;0.999;1]
%for  uuu=1:1:length(RR)
disp('RI')
tic
Rx=0.99;
Ry=0.99;
Rz=0.99;
%ordre=100 %ordre de calcul... nombre de r?flexions envisag?es
N=60000
d=1; %dirac
for ee=1:MM
    X_1=l*rand;
    Y_1=p*rand;
    Z_1=h*rand;
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
%                           Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
%                                     -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
%                                     -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];
%             %
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
            
%                         Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
%                                     (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
%                                     (cb*sa) (sb*sa) ca];
            
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
    SSz(:,ee)=Az;
    SSx(:,ee)=Ax;
    SSy(:,ee)=Ay;
    
end

filename5 = sprintf('resultietr_angle%d.mat',tilt);

save(filename5,'freq','SSx','SSy','SSz')
clear
% figure(1)
% % hold on
% plot(freq,cumsum(RRfinalx),freq,cumsum(RRfinaly),freq,cumsum(RRfinalz))%,freq,m2*freq+p2,'k',freq,m1*freq+p1,'k')
% ylabel('rejets cumul???s')
% xlabel('Fr???quence')
% title([num2str(sum(RRfinalx+RRfinaly+RRfinalz)),' rejets pour ',num2str(3*length(freq)),' tests, ',num2str(MM),' ???chantillons,','R=',num2str(Rx),', ',num2str(l),'x',num2str(p),'x',num2str(h)])
% grid on
% xlim([0 max(freq)])
