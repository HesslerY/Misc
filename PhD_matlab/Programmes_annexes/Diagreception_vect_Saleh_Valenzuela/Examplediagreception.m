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
Lt = 6e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9)));
load(filename)
toc
%Loss coefficient
R = .998;

%Reception point rectangular coordinates
X_1 = 6.5;
Y_1 = 3;
Z_1 = 1.5;

[Sx,Sy,Sz,t,azim,elev] = CIRvect(X_1,Y_1,Z_1);
toc
D=2% taille du syst�me en lambda

n = 5;
xa = D*rand(n,1);
ya = D*rand(n,1);
za = D*rand(n,1);
Ia = complex(rand(n,1),rand(n,1));
S = 2;
fr = fres(elev, azim, n, xa, ya, za, Ia);
fa = fantres(fr, elev, azim, S);


Signalx=Sx.*fa;
Signaly=Sy.*fa;
Signalz=Sz.*fa;


Snx=zeros(1,N);
Sny=zeros(1,N);
Snz=zeros(1,N);
Snxa=zeros(1,N);
Snya=zeros(1,N);
Snza=zeros(1,N);
zn=round((N-1)*t./Lt);

for j=1:1:length(Signalx)
    if zn(j)<N+1
        Snx(zn(j))=Snx(zn(j))+Sx(j);
        Sny(zn(j))=Sny(zn(j))+Sy(j);
        Snz(zn(j))=Snz(zn(j))+Sz(j);
        
        Snxa(zn(j))=Snxa(zn(j))+Signalx(j);
        Snya(zn(j))=Snya(zn(j))+Signaly(j);
        Snza(zn(j))=Snza(zn(j))+Signalz(j);
    end
end
clear Signalx Signaly Signalz
toc

figure(1)
plot(tt,Snz,tt,Snza./sqrt(sum(Snza.^2)/sum(Snz.^2)))

%Pulsed signal
tau = .3e-6; %length of the pulse in seconds
f0 = 1e9; %monochromatic pulse frequency

N = round(10*Lt*f0) %number of points for the chosen time-window (Lt)
tt = 0:Lt/(N-1):Lt; %time scale
x = 0:1/((N-1)/Lt):tau;
s = sin(2*pi*f0*x); %pulsed signal


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