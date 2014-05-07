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
Lt = 1e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)

%Pulsed signal
tau = Lt; %length of the pulse in seconds
f0 = 1e9; %monochromatic pulse frequency

l=8.7;
p=3.7;
h=2.9;

N = round(10*Lt*f0) %number of points for the chosen time-window (Lt)
t = 0:Lt/(N-1):Lt; %time scale
x = 0:1/((N-1)/Lt):tau;
s = sin(2*pi*f0*x); %pulsed signal
X=0:.01:l
%Loss coefficient 
R = 0.998;
for i=1:length(X);
    disp(i)
%Reception point rectangular coordinates
X_1 = X(i);
Y_1 = 2;
Z_1 = 1.5;

[Sx,Sy,Sz] = CIR(X_1,Y_1,Z_1);
ssz=conv(s,Sz);

signal(i,:)=ssz(1:N);
end



imagesc(X,t/1e-9,10*log10(signal'.^2))
caxis([-30 30])
%saveas(gcf,'TPsquelettedB.eps','epsc2')
