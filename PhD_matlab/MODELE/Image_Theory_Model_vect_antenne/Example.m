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

global Lt c Rx Ry Rz N POS

tic
freq=[];

c=3e8;%
Lt=1e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance
%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns.mat',nbre_elements,round(Lt/(1e-9))); 
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
[Sx,Sy,Sz]=CIR(X_1,Y_1,Z_1);
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
