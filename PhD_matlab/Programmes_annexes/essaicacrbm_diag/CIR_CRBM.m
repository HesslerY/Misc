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

    CONVomni1=[];
    CONVant1=[];
    CONVomni5=[];
    CONVant5=[];
    CONVomni10=[];
    CONVant10=[];
    CONVomni10dB=[];
    CONVant10dB=[];
    

load XYZ150b.mat
load angles.mat
for uuu=1:150 
    disp(uuu/150)
    

    
   
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    load antenneD1n5.mat
    %  load antenneD5n5.mat
    %  load antenneD10n5.mat
    %  load antenne10dB.mat
    
    
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
    
    CONVomni1=[CONVomni1;sigomni(1:N)];
    CONVant1=[CONVant1;sigantenne(1:N)];

    save('essaiCRBMD1n5.mat','tt','CONVomni1','CONVant1');
toc
    

    
   
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    %load antenneD1n5.mat
      load antenneD5n5.mat
    %  load antenneD10n5.mat
    %  load antenne10dB.mat
    
    
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
    
    CONVomni5=[CONVomni5;sigomni(1:N)];
    CONVant5=[CONVant5;sigantenne(1:N)];

    save('essaiCRBMD5n5.mat','tt','CONVomni5','CONVant5');
    toc
    

    
   
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    %  load antenneD1n5.mat
    %  load antenneD5n5.mat
       load antenneD10n5.mat
    %  load antenne10dB.mat
    
    
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
    
    CONVomni10=[CONVomni10;sigomni(1:N)];
    CONVant10=[CONVant10;sigantenne(1:N)];

    save('essaiCRBMD10n5.mat','tt','CONVomni10','CONVant10');
    
    toc
    

   
    filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
    load(filename)
    %load antenneD1n5.mat
    %  load antenneD5n5.mat
    %  load antenneD10n5.mat
      load antenne10dB.mat
    
    
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
    
    CONVomni10dB=[CONVomni10dB;sigomni(1:N)];
    CONVant10dB=[CONVant10dB;sigantenne(1:N)];

    save('essaiCRBMD10dB.mat','tt','CONVomni10dB','CONVant10dB');
    toc
end


