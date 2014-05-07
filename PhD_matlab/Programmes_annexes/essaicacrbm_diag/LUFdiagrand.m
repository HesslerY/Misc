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
for dx=1:5
FFTx=[];
FFTy=[];
FFTz=[];


FFTxa=[];
FFTya=[];
FFTza=[];

c = 299792458;%
Lt = 6e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8e111.mat',nbre_elements,round(Lt/(1e-9)));
load(filename)
toc
%Loss coefficient
R = .998;
tau=-3/(2*c*log(R))
%Reception point rectangular coordinates

    load XYZ150b.mat
    for uuu=1:50
        disp(uuu)
        [Sx,Sy,Sz,t,azim,elev] = CIRvect(X(uuu),Y(uuu),Z(uuu));
        D=5;% taille du syst�me en lambda
        
        
n =2;

xa = 1:dx:dx*n+1;%D*rand(n,1);
ya = ones(1,n);%D*rand(n,1);
za = ones(1,n);%D*rand(n,1);

%[Ia, dph] = dolph3(1, dx, n,10)
Ia=ones(1,n)
dph=zeros(1,n-1)
S = 2;
        
        %Pulsed signal
        tau = Lt; %length of the pulse in seconds
        f0 = 1e9; %monochromatic pulse frequency
        
        N = round(10*Lt*f0) %number of points for the chosen time-window (Lt)
        tt = 0:Lt/(N-1):Lt; %time scale
        x = 0:1/((N-1)/Lt):tau;
        s = sin(2*pi*f0*x); %pulsed signal
        
        fr = fres(mod(pi/2-elev+pi/2-acos(sqrt(2/3)),pi), mod(pi-azim+pi/4,2*pi), n, xa, ya, za, Ia);
        fa = fantres(fr,mod(pi/2-elev+pi/2-acos(sqrt(2/3)),pi), mod(pi-azim+pi/4,2*pi), S);
        
        clear fr
        
        clear azim elev
        
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
        
        clear Sx Sy Sz Signalx Signaly Signalz t
        
        %sigconvz=conv(s,Snza);
        %CIRz=[CIRz;Snza];
        %CONVz=[CONVz;sigconvz(1:N)];
        
        Fs = N/Lt; %sampling frequency
        T = 1/Fs;                     % sample length
        L = N;                  %Number of points
        tt = (0:L-1)*T;                % time...
        
        NFFT = 2^nextpow2(N); % Next power of 2 from length of y
        Yx = fft((Snx),NFFT);
        Yy = fft((Sny),NFFT);
        Yz = fft((Snz),NFFT);
        
        Yxa = fft((Snxa),NFFT);
        Yya = fft((Snya),NFFT);
        Yza = fft((Snza),NFFT);
        f = Fs/2*linspace(0,1,NFFT/2);
        
        % Three axis frequency response
        FFTx(uuu,:)=Yx(1:NFFT/2);
        FFTy(uuu,:)=Yy(1:NFFT/2);
        FFTz(uuu,:)=Yz(1:NFFT/2);
        
        FFTxa(uuu,:)=Yxa(1:NFFT/2);
        FFTya(uuu,:)=Yya(1:NFFT/2);
        FFTza(uuu,:)=Yza(1:NFFT/2);
        freq=f;
        %if mod(uuu,50)==0
        %end
        toc
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Diag
nt = 100;
np = 191;
theta2 = linspace(0, pi, nt);
phi2 = linspace(0, 2*pi, np);
[thetaa, phia] = meshgrid(theta2, phi2);
frant = fres(thetaa, phia, n, xa, ya, za, Ia);
faant = fantres(frant, thetaa, phia, S);
d1ant = Directivite(faant, thetaa, phia);
%tracediag(faant,thetaa, phia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D=d1ant; 
    filenamee=sprintf('FFT_50-6muantenneorientationdx_%d_elem.mat',dx)
    save(filenamee,'N','FFTx','FFTy','FFTz','f','FFTxa','FFTya','FFTza','D')
end