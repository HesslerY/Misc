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
% global Lt c R POS
% 
% % FFTx=[];
% % FFTy=[];
% % FFTz=[];
% % 
% % FFTxa=[];
% % FFTya=[];
% % FFTza=[];

c = 299792458;%
Lt = 6e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance
CONVomni=[];
CONVant=[];
%loading the Position matrix from the image generator
%filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
%load(filename)

% POS=[1 2 1 0 0 0];
% 
% toc
% %Loss coefficient
% R = .995;
% tau=-3/(2*c*log(R))
% %Reception point rectangular coordinates
% 
% 
% % load antenneD5n5.mat
% % load antenneD10n5.mat
% % load antenne10dB.mat
% 
% %    load XYZ150b.mat
% 
% %    disp([num2str(ang),'/50'])
% 
% 
% [Sx,Sy,Sz,t,azim,elev] = CIRvectca(1,5,1);
% 
% save('CACIR6mu.mat','Sx','Sy','Sz','azim','elev','t')

load 'CACIR6mu.mat'

load antenne10dB.mat

alp=pi*rand(1,150);
bet=2*pi*rand(1,150);

%Pulsed signal
tau = Lt; %length of the pulse in seconds
f0 = 1e9; %monochromatic pulse frequency

N = round(10*Lt*f0) %number of points for the chosen time-window (Lt)
tt = 0:Lt/(N-1):Lt; %time scale
x = 0:1/((N-1)/Lt):tau;
s = sin(2*pi*f0*x); %pulsed signal

 fr = fres(mod(alp,pi), mod(bet,2*pi), n, xa, ya, za, Ia);
 fa = fantres(fr,mod(alp,pi), mod(bet,2*pi), S);
 
 clear fr
 clear azim elev
  
 for ang=1:150
    Signalx=Sx.*fa(ang);
    Signaly=Sy.*fa(ang);
    Signalz=Sz.*fa(ang);
        
    Snx=zeros(1,N);
    Sny=zeros(1,N);
    Snz=zeros(1,N);
    Snxa=zeros(1,N);
    Snya=zeros(1,N);
    Snza=zeros(1,N);
    zn=round((N-1)*t./Lt);
    
    for j=1:1:length(Signalz)
        if zn(j)<N+1
            Snx(zn(j))=Snx(zn(j))+Sx(j);
            Sny(zn(j))=Sny(zn(j))+Sy(j);
            Snz(zn(j))=Snz(zn(j))+Sz(j);
            
            Snxa(zn(j))=Snxa(zn(j))+Signalx(j);
            Snya(zn(j))=Snya(zn(j))+Signaly(j);
            Snza(zn(j))=Snza(zn(j))+Signalz(j);
        end
    end
    Sigomni=Snx*cos(bet(ang))*sin(alp(ang))+Sny*sin(bet(ang))*sin(alp(ang))+Snz*cos(alp(ang));
    Siga=Snxa*cos(bet(ang))*sin(alp(ang))+Snya*sin(bet(ang))*sin(alp(ang))+Snza*cos(alp(ang));
    %Sigomni=Snz;
    %Siga=Snza;
    %clear Sx Sy Sz Signalx Signaly Signalz t
    sigomni=conv(s,Sigomni);
    sigantenne=conv(s,Siga);
    %CIRz=[CIRz;Snza];
    CONVomni=[CONVomni;sigomni(1:N)];
    CONVant=[CONVant;sigantenne(1:N)];

    toc
end
save('essaiCAD10dB.mat','tt','CONVomni','CONVant','alp','bet');