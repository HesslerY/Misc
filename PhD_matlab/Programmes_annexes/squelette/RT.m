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
Lt = 12e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8RT.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)
toc
%Pulsed signal
tau = 1e-6; %length of the pulse in seconds
f0 = 1e9; %monochromatic pulse frequency

N = round(5*Lt*f0) %number of points for the chosen time-window (Lt)
t = 0:Lt/(N-1):Lt; %time scale
x = 0:1/((N-1)/Lt):tau;
s = sin(2*pi*f0*x); %pulsed signal

%Loss coefficient 
R = 0.998;

%Reception point rectangular coordinates
X_1 = 4.5;
Y_1 = 3;
Z_1 = 1.5;

[Sx1,Sy1,Sz1] = CIR(X_1,Y_1,Z_1);
toc
%Reception point rectangular coordinates
X_2 = 4.1;
Y_2 = 3;
Z_2 = 1.5;

[Sx2,Sy2,Sz2] = CIR(X_2,Y_2,Z_2);
toc

%Reception point rectangular coordinates
X_3 = 4.25;
Y_3 = 3;
Z_3 = 1.5;

[Sx3,Sy3,Sz3] = CIR(X_3,Y_3,Z_3);
toc

save('RT.mat','t','Sx1','Sy1','Sz1','Sx2','Sy2','Sz2','Sx3','Sy3','Sz3')





