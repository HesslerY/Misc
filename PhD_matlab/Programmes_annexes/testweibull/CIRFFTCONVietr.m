clear

tic
freq=[];

c=3e8;%
Lt=3e-6;
%
%ordre=201;%round(dmax/min
dmax=c*Lt
filename = sprintf('1elemIETR_3000ns.mat');

load(filename)
N=1000000
%Signal puls�
tau=.3e-6; %longueur du pulse en s
f0=1e9; %µfréquence du pulse en Hz
%N=round(30*Lt*f0) %nombre de points de la simulation globale

t=0:Lt/(N-1):Lt; %�chelle de temps

x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %train d'onde sinusoidal
%N=round(30*Lt*f0) %nombre de points de la simulation globale

length(POS);
%RR=[0.8;0.90;0.92;0.95;0.98;0.99;0.995;0.999;1]
%for  uuu=1:1:length(RR)
disp('RI')
tic
Rx=0.998;
Ry=0.998;
Rz=0.998;
%ordre=100 %ordre de calcul... nombre de r?flexions envisag?es

d=1; %dirac

    X_1=6;
    Y_1=2;
    Z_1=1.5;
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
    end    
    toc
    
%Convolutions
disp('CONV')
Six=conv(Sx,s);    
Siy=conv(Sy,s);
Siz=conv(Sz,s);    
  Signalfinalx=Six(1:N);
  Signalfinaly=Siy(1:N);
  Signalfinalz=Siz(1:N);  
toc
    disp('FFT')
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
    freq=f;
    FFTz=Az;
    FFTx=Ax;
    FFTy=Ay;
   toc 


filename5 = sprintf('RIFFTCONV.mat');

save(filename5,'t','freq','Sx','Sy','Sz','Signalfinalx','Signalfinaly','Signalfinalz','FFTx','FFTy','FFTz')
