%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                               %%%%%%
%%%%%       CHANNEL IMPULSE RESPONSE, PULSED SIGNAL RESPONSE        %%%%%%
%%%%%                    & FREQUENCY RESPONSE                       %%%%%%
%%%%%              by E. Amador (manuamador@gmail.com)              %%%%%%
%%%%%                         IETR/DGA                              %%%%%%
%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear


tic
freq=[];

c=3e8;%
Lt=1e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance
%loading the Position matrix from the image generator
filename = sprintf('%delem_%dnslongue.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)

%Pulsed signal
tau=.3e-6; %length of the pulse in seconds
f0=1e9; %monochromatic pulse frequency

N=round(30*Lt*f0) %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.998;
Ry=0.998;
Rz=0.998;


%Reception point coordinates
X_1=4.5;
Y_1=3;
Z_1=1.5;


disp('CIR')
d=1; %value of the elementary pulse
tic

%Channel Impulse Response
Sx=zeros(1,N);
Sy=zeros(1,N);
Sz=zeros(1,N);
%for j=1:1:length(POS) %loop over the image-currents... to be vctorized... one day may be...
    
    DX=X_1-POS(:,1);
    DY=Y_1-POS(:,2);
    DZ=Z_1-POS(:,3);
    
    dist=sqrt(DX.^2+DY.^2+DZ.^2);
    zl=dist./c;
    zn=round((N-1)*dist./c/Lt);
    
    
        alpha=POS(:,9);
        beta=POS(:,10);
        ca=cos(alpha);
        sa=sin(alpha);
        cb=cos(beta);
        sb=sin(beta);
        
        %Ralpha/beta, coordinates transformation matrix:
        %                           Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
        %                                     -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
        %                                     -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];
        
        %rectangular coordinates calculation in the local system attached to the considered current (developped expressions, a matrix product is slower)
        distx=((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
        disty=(-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
        distz=(-cb.*sa).*DX+(-sb.*sa).*DY+(ca).*DZ;
        
        distxy=sqrt(distx.^2+disty.^2);
        
       
        E=Rx.^POS(:,4).*Ry.^POS(:,5).*Rz.^POS(:,6).*POS(:,8).*d./dist; % Free-space attenuation
        costheta=distz./dist;
        sintheta=distxy./dist;
        cosphi=distx./distxy;
        sinphi=disty./distxy;
        Antth=-sintheta; %dipole radiation pattern
        
        %Reverse transformation matrix
        %                         Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
        %                                     (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
        %                                     (cb*sa) (sb*sa) ca];
        
      
    
        
        %Three axis channel impulse response construction
        Sxb=E.*(((-sb).^2+(1-(-sb).^2).*ca).*(Antth.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(Antth.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*Antth));
        Syb=E.*((-sb.*cb.*(1-ca)).*(Antth.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(Antth.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*Antth));
        Szb=E.*((cb.*sa).*(Antth.*costheta.*cosphi)+(sb.*sa).*(Antth.*costheta.*sinphi)+ca.*(-sintheta.*Antth));

toc
Sx=zeros(1,N);
Sy=zeros(1,N);
Sz=zeros(1,N);


for j=1:1:length(Sxb)
    if zn(j)<N+1
Sx(zn(j))=Sx(zn(j))+Sxb(j);
Sy(zn(j))=Sy(zn(j))+Syb(j);
Sz(zn(j))=Sz(zn(j))+Szb(j);
    end
end
toc

%Convolution of the CIRs with the chosen pulsed signal
disp('CONV')
Six=conv(Sx,s);
Siy=conv(Sy,s);
Siz=conv(Sz,s);
Signalfinalx=Six(1:N);
Signalfinaly=Siy(1:N);
Signalfinalz=Siz(1:N);
toc

%Frequency response
disp('FFT')


Fs = N/Lt; %sampling frequency
T = 1/Fs;                     % sample length
L = N;                  %Number of points
tt = (0:L-1)*T;                % time...

NFFT = 2^nextpow2(N); % Next power of 2 from length of y
Yx = fft((Sx),NFFT)/N;
Yy = fft((Sy),NFFT)/N;
Yz = fft((Sz),NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2);

% Three axis frequency response
FFTx=abs(Yx(1:NFFT/2));
FFTy=abs(Yy(1:NFFT/2));
FFTz=abs(Yz(1:NFFT/2));
freq=f;

toc

%CIR figure
figure(1)
subplot(3,1,1)
plot(t,Sx)
title('E_x')
grid on
xlabel('time in s')
ylabel('V/m')

subplot(3,1,2)
plot(t,Sy)
title('E_y')
grid on
xlabel('time in s')
ylabel('V/m')

subplot(3,1,3)
plot(t,Sz)
title('E_z')
grid on
xlabel('time in s')
ylabel('V/m')

%Pulsed signal figure
figure(2)
subplot(3,1,1)
plot(t,Signalfinalx)
title('E_x')
grid on
xlabel('time in s')
ylabel('V/m')

subplot(3,1,2)
plot(t,Signalfinaly)
title('E_y')
grid on
xlabel('time in s')
ylabel('V/m')

subplot(3,1,3)
plot(t,Signalfinalz)
title('E_z')
grid on
xlabel('time in s')
ylabel('V/m')

%Frequency response figure
figure(3)
subplot(3,1,1)
plot(freq/1e6,20*log10(FFTx))
title('FFT_x')
xlim([0 500])
grid on
xlabel('frequency in MHz')
ylabel('dB')

subplot(3,1,2)
plot(freq/1e6,20*log10(FFTy))
title('FFT_y')
xlim([0 500])
grid on
xlabel('frequency in MHz')
ylabel('dB')

subplot(3,1,3)
plot(freq/1e6,20*log10(FFTz))
title('FFT_z')
xlim([0 500])
grid on
xlabel('frequency in MHz')
ylabel('dB')
