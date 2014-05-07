%Example.m sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       CHANNEL IMPULSE RESPONSE, PULSED SIGNAL RESPONSE        %
%                    & FREQUENCY RESPONSE                       %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
tic
global Lt c R POS


c = 299792458;%
Lt = .04e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8e.mat',nbre_elements,round(Lt/(1e-9)));
load(filename)
toc
%Loss coefficient
R = .998;

%Reception point rectangular coordinates
X_1 = 5.28;
Y_1 = 2.29;
Z_1 = 1.01;

[Sx,Sy,Sz,t,azim,elev] = CIRvect(X_1,Y_1,Z_1);
toc

% D=5% taille du syst�me en lambda
% 
% n = 50;
% xa = D*rand(n,1);
% ya = D*rand(n,1);
% za = D*rand(n,1);
% Ia = complex(rand(n,1),rand(n,1));
% S = 2;
% save('ant.mat','xa','ya','za','Ia','S','n')
% 
% load('ant.mat')

n =3;
dx=0.25
xa = 1:dx:dx*n+1;%D*rand(n,1);
ya = ones(1,n);%D*rand(n,1);
za = ones(1,n);%D*rand(n,1);

[Ia, dph] = dolph3(1, dx, n, 10)


S = 2;

%Pulsed signal
tau = Lt; %length of the pulse in seconds
f0 = 1e9; %monochromatic pulse frequency

N = round(10*Lt*f0) %number of points for the chosen time-window (Lt)
tt = 0:Lt/(N-1):Lt; %time scale
x = -tau/2:1/((N-1)/Lt):tau/2;
%s = sin(2*pi*f0*x); %pulsed signal
s=gauspuls(x,882e6,0.3)

Snz=zeros(1,N);
AZIM=zeros(1,N);
ELEV=zeros(1,N);
zn=round((N-1)*t./Lt);
for j=1:1:length(Sz)
    if zn(j)<N+1
%         Snx(zn(j))=Snx(zn(j))+Sx(j);
%         Sny(zn(j))=Sny(zn(j))+Sy(j);
        Snz(zn(j))=Snz(zn(j))+Sz(j);
        ELEV(zn(j))=-elev(j)+pi/2;
        AZIM(zn(j))=-azim(j)+pi;
    end
end

plot(x,s)

CIRz=[];
CONVz=[];
FFTz=[];
M=360;
for o=2*pi/M:2*pi/M:2*pi
disp(o/(2*pi))
azimp=azim+o;    
fr = fres(elev, azimp, n, xa, ya, za, Ia);
fa = fantres(fr, elev, azimp, S);


% Signalx=Sx.*fa;
% Signaly=Sy.*fa;
Signalz=Sz.*fa;



% Snxa=zeros(1,N);
% Snya=zeros(1,N);
Snza=zeros(1,N);

for j=1:1:length(Signalz)
    if zn(j)<N+1
%         Snxa(zn(j))=Snxa(zn(j))+Signalx(j);
%         Snya(zn(j))=Snya(zn(j))+Signaly(j);
        Snza(zn(j))=Snza(zn(j))+Signalz(j);

    end
end
sigconvz=conv(Snza,s);

CIRz=[CIRz;Snza];
CONVz=[CONVz;sigconvz(1:N)];

Fs = N/Lt; %sampling frequency
T = 1/Fs;                     % sample length
L = N;                  %Number of points
tt = (0:L-1)*T;                % time...

NFFT = 2^nextpow2(N); % Next power of 2 from length of y
%Yx = fft((Sx),NFFT)/N;
%Yy = fft((Sy),NFFT)/N;
Yz = fft((Snza),NFFT);
f = Fs/2*linspace(0,1,NFFT/2);

% Three axis frequency response
%FFTx=abs(Yx(1:NFFT/2));
%FFTy=abs(Yy(1:NFFT/2));
FFTz=[FFTz;abs(Yz(1:NFFT/2))];
freq=f;

clear Signalx Signaly Signalz
toc
end

%for i=1:360
for i=[110,180,40,250,360]
%figure(1)
% plot(tt,abs(Snz/max(Snz)),tt,abs(CIRz(i,:)/max(abs(CIRz(i,:)))), 'DisplayName', 'Snza', 'YDataSource', 'Snza'); 
% title(num2str(i))
% xlim([10e-9 20e-9])
% ylim([0 1.1])
% xlabel('time [s]')
% ylabel('[V/m]')
% grid on

U=find(abs(ELEV-pi/2)<pi/10);

elevh=ELEV(U);
azimh=azim(U);
CIRh140=CIRz(i,U);%/max(abs(CIRz(140,U)));
figure,
plot(azimh*180/pi,CIRh140.^2,'.')
title(num2str(i))

getframe;

end


U=find(abs(ELEV-pi/2)<pi/5);

elevh=ELEV(U);
azimh=azim(U);
CIRh140=CIRz(180,U);%/max(abs(CIRz(140,U)));

plot(azim*180/pi,CIRh140.^2)




% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Diag
nt = 100;
np = 191;
theta2 = linspace(0, pi, nt);
phi2 = linspace(0, 2*pi, np);
[thetaa, phia] = meshgrid(theta2, phi2);
frant = fres(thetaa, phia, n, xa, ya, za, Ia);
faant = fantres(frant, thetaa, phia, S);
d1ant = Directivite(faant, thetaa, phia)
tracediag(faant,thetaa, phia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%