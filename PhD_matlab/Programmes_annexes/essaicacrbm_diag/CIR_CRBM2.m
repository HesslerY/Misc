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
CONVomni=[];
CONVant=[];
%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
load(filename)

%POS=[1 2 1 0 0 0];


%Loss coefficient
R = .995;

CONVomni7dB=[];
CONVant7dB=[];
CONVomni17dB=[];
CONVant17dB=[];
CONVomni21dB=[];
CONVant21dB=[];


load XYZ150b.mat
load angles.mat
for uuu=1:150
    disp(uuu/150)
    
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    
    load antenne7dB.mat
    
    
    [Sx,Sy,Sz,t,azim,elev] = CIRvect(X(uuu),Y(uuu),Z(uuu));
    %filename=sprintf('CRBMCIR6mu_%d.mat',uuu)
    %save(filename,'Sx','Sy','Sz','azim','elev','t')
    clear POS
    
    %Pulsed signal
    tau = Lt; %length of the pulse in seconds
    f0 = 1e9; %monochromatic pulse frequency
    
    N = round(10*Lt*f0); %number of points for the chosen time-window (Lt)
    tt = 0:Lt/(N-1):Lt; %time scale
    x = 0:1/((N-1)/Lt):tau;
    s = sin(2*pi*f0*x); %pulsed signal
    
    fr = fres(mod(pi/2-elev-alp(uuu),pi), mod(pi-azim-bet(uuu),2*pi), n, xa, ya, za, Ia);
    fa = fantres(fr,mod(pi/2-elev-alp(uuu),pi),mod(pi-azim-bet(uuu),2*pi), S);
    clear fr azim elev
    
    Signalx=Sx.*fa;
    Signaly=Sy.*fa;
    Signalz=Sz.*fa;
    
    clear fa
    
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
    
    clear Sx Sy Sz Signalx Signaly Signalz
    
    Sigomni=Snx*cos(bet(uuu))*sin(alp(uuu))+Sny*sin(bet(uuu))*sin(alp(uuu))+Snz*cos(alp(uuu));
    Siga=Snxa*cos(bet(uuu))*sin(alp(uuu))+Snya*sin(bet(uuu))*sin(alp(uuu))+Snza*cos(alp(uuu));
    
    sigomni=conv(s,Sigomni);
    sigantenne=conv(s,Siga);
    
    CONVomni7dB=[CONVomni7dB;sigomni(1:N)];
    CONVant7dB=[CONVant7dB;sigantenne(1:N)];
    
    save('essaiCRBMD7dB.mat','tt','CONVomni7dB','CONVant7dB');
    
    
    
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    %load antenneD1n5.mat
    %  load antenneD5n5.mat
    %  load antenneD10n5.mat
    load antenne17dB.mat
    
    
    [Sx,Sy,Sz,t,azim,elev] = CIRvect(X(uuu),Y(uuu),Z(uuu));
    %filename=sprintf('CRBMCIR6mu_%d.mat',uuu)
    %save(filename,'Sx','Sy','Sz','azim','elev','t')
    clear POS
    
    %Pulsed signal
    tau = Lt; %length of the pulse in seconds
    f0 = 1e9; %monochromatic pulse frequency
    
    N = round(10*Lt*f0); %number of points for the chosen time-window (Lt)
    tt = 0:Lt/(N-1):Lt; %time scale
    x = 0:1/((N-1)/Lt):tau;
    s = sin(2*pi*f0*x); %pulsed signal
    
    fr = fres(mod(pi/2-elev-alp(uuu),pi), mod(pi-azim-bet(uuu),2*pi), n, xa, ya, za, Ia);
    fa = fantres(fr,mod(pi/2-elev-alp(uuu),pi),mod(pi-azim-bet(uuu),2*pi), S);
    clear fr azim elev
    
    Signalx=Sx.*fa;
    Signaly=Sy.*fa;
    Signalz=Sz.*fa;
    
    clear fa
    
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
    
    clear Sx Sy Sz Signalx Signaly Signalz
    
    Sigomni=Snx*cos(bet(uuu))*sin(alp(uuu))+Sny*sin(bet(uuu))*sin(alp(uuu))+Snz*cos(alp(uuu));
    Siga=Snxa*cos(bet(uuu))*sin(alp(uuu))+Snya*sin(bet(uuu))*sin(alp(uuu))+Snza*cos(alp(uuu));
    
    sigomni=conv(s,Sigomni);
    sigantenne=conv(s,Siga);
    
    CONVomni17dB=[CONVomni17dB;sigomni(1:N)];
    CONVant17dB=[CONVant17dB;sigantenne(1:N)];
    
    save('essaiCRBMD17dB.mat','tt','CONVomni17dB','CONVant17dB');
    
    
    toc
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    %load antenneD1n5.mat
    %  load antenneD5n5.mat
    %  load antenneD10n5.mat
    load antenne21dB.mat
    
    
    [Sx,Sy,Sz,t,azim,elev] = CIRvect(X(uuu),Y(uuu),Z(uuu));
    %filename=sprintf('CRBMCIR6mu_%d.mat',uuu)
    %save(filename,'Sx','Sy','Sz','azim','elev','t')
    clear POS
    
    %Pulsed signal
    tau = Lt; %length of the pulse in seconds
    f0 = 1e9; %monochromatic pulse frequency
    
    N = round(10*Lt*f0); %number of points for the chosen time-window (Lt)
    tt = 0:Lt/(N-1):Lt; %time scale
    x = 0:1/((N-1)/Lt):tau;
    s = sin(2*pi*f0*x); %pulsed signal
    
    fr = fres(mod(pi/2-elev-alp(uuu),pi), mod(pi-azim-bet(uuu),2*pi), n, xa, ya, za, Ia);
    fa = fantres(fr,mod(pi/2-elev-alp(uuu),pi),mod(pi-azim-bet(uuu),2*pi), S);
    clear fr azim elev
    
    Signalx=Sx.*fa;
    Signaly=Sy.*fa;
    Signalz=Sz.*fa;
    
    clear fa
    
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
    
    clear Sx Sy Sz Signalx Signaly Signalz
    
    Sigomni=Snx*cos(bet(uuu))*sin(alp(uuu))+Sny*sin(bet(uuu))*sin(alp(uuu))+Snz*cos(alp(uuu));
    Siga=Snxa*cos(bet(uuu))*sin(alp(uuu))+Snya*sin(bet(uuu))*sin(alp(uuu))+Snza*cos(alp(uuu));
    
    sigomni=conv(s,Sigomni);
    sigantenne=conv(s,Siga);
    
    CONVomni21dB=[CONVomni21dB;sigomni(1:N)];
    CONVant21dB=[CONVant21dB;sigantenne(1:N)];
    
    save('essaiCRBMD21dB.mat','tt','CONVomni21dB','CONVant21dB');
    toc
end


