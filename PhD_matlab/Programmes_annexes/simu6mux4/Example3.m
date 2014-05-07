%Example.m sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       CHANNEL IMPULSE RESPONSE, PULSED SIGNAL RESPONSE        %
%                    & FREQUENCY RESPONSE                       %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

global Lt c R N POS

tic
c = 299792458;%
Lt = 6e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)
l=8.7;
p=3.7;
h=2.9;
%Pulsed signal
tau = 1e-6; %length of the pulse in seconds
f0 = 2e9; %monochromatic pulse frequency

N = round(5*Lt*f0) %number of points for the chosen time-window (Lt)
t = 0:Lt/(N-1):Lt; %time scale
x = 0:1/((N-1)/Lt):tau;
s = sin(2*pi*f0*x); %pulsed signal

%Loss coefficient 
R = 0.998;

%NN=50 %number of positions

%lambda0=1.5; %wavelength
%1/8
%first position
%X=[(l-lambda0)*rand+lambda0/2];
%Y=[(p-lambda0)*rand+lambda0/2];
%Z=[(h-lambda0)*rand+lambda0/2];

%while (length(X)<NN) %the other positions
 %   X_1=(l-lambda0)*rand+lambda0/2;
  %  Y_1=(p-lambda0)*rand+lambda0/2;
  %  Z_1=(h-lambda0)*rand+lambda0/2;
  %  D=sqrt((X_1-X).^2+(Y_1-Y).^2+(Z_1-Z).^2);
  %  if min(D)>lambda0/2
  %      X=[X;X_1];
  %      Y=[Y;Y_1];
  %      Z=[Z;Z_1];
  %  end
%end
load XYZ150b.mat

disp('OK')
%Reception point rectangular coordinates
% X_1 = 4.5;
% Y_1 = 3;
% Z_1 = 1.5;
for uuu=76:113
    X_1=X(uuu);
    Y_1=Y(uuu);
    Z_1=Z(uuu);
 tic   
[Sx,Sy,Sz] = CIR(X_1,Y_1,Z_1);
disp(uuu/400)
% %Convolution of the CIRs with the chosen pulsed signal
% disp('Convolution...')
 Six = conv(Sx,s);
% Siy = conv(Sy,s);
 Siz = conv(Sz,s);
 Signalfinalx(uuu,:) = Six(1:N);
% Signalfinaly = Siy(1:N);
 Signalfinalz(uuu,:) = Siz(1:N);
% toc

%Frequency response
%disp('FFT...')
%Fs = N/Lt; %sampling frequency
%T = 1/Fs;                     % sample length
%L = N;                  %Number of points
%tt = (0:L-1)*T;                % time...

%NFFT = 2^nextpow2(N); % Next power of 2 from length of y
%Yx = fft((Sx),NFFT)/N;
%Yy = fft((Sy),NFFT)/N;
%Yz = fft((Sz),NFFT)/N;
%f = Fs/2*linspace(0,1,NFFT/2);

% Three axis frequency response
%FFTx(uuu,:) = (Yx(1:NFFT/2));
%FFTy(uuu,:) = (Yy(1:NFFT/2));
%FFTz(uuu,:) = (Yz(1:NFFT/2));
%freq = f;
%toc
end
save('result3_6MU.mat','Signalfinalx','Signalfinalz','t','R','tau')
% 
% %CIR figure
% figure(1)
% subplot(3,1,1)
% plot(t/1e-6,Sx)
% title('E_x')
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% ylim([-0.2 0.2])
% subplot(3,1,2)
% plot(t/1e-6,Sy)
% title('E_y')
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% ylim([-0.2 0.2])
% subplot(3,1,3)
% plot(t/1e-6,Sz)
% title('E_z')
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% ylim([-0.2 0.2])
% %Pulsed signal figure
% figure(2)
% subplot(3,1,1)
% plot(t/1e-6,Signalfinalx)
% title('E_x')
% ylim([-6 +6])
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% 
% subplot(3,1,2)
% plot(t/1e-6,Signalfinaly)
% title('E_y')
% ylim([-6 +6])
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% 
% subplot(3,1,3)
% plot(t/1e-6,Signalfinalz)
% title('E_z')
% ylim([-6 +6])
% grid on
% xlabel('temps en $\mu$s')
% ylabel('V/m')
% 
% %Frequency response figure
% figure(3)
% subplot(3,1,1)
% plot(freq/1e6,20*log10(FFTx))
% title('FFT_x')
% xlim([0 500])
% ylim([-150 -50])
% grid on
% xlabel('Frequence en MHz')
% ylabel('dB/Hz')
% 
% subplot(3,1,2)
% plot(freq/1e6,20*log10(FFTy))
% title('FFT_y')
% xlim([0 500])
% ylim([-150 -50])
% grid on
% xlabel('Frequence en MHz')
% ylabel('dB/Hz')
% 
% subplot(3,1,3)
% plot(freq/1e6,20*log10(FFTz))
% title('FFT_z')
% xlim([0 500])
% ylim([-150 -50])
% grid on
% xlabel('Frequence en MHz')
% ylabel('dB/Hz')
