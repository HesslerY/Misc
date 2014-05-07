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

freq=[];
l=8.7;
p=3.7;
h=2.9;

c=3e8;%
Lt=1e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance

tau=.3e-6; %length of the pulse in seconds
f0=.5e9; %monochromatic pulse frequency

N=600000 %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.998;
Ry=0.998;
Rz=0.998;


%Reception point coordinates

lambda=c/f0;

X=[(l-lambda)*rand+lambda/2];
Y=[(p-lambda)*rand+lambda/2];
Z=[(h-lambda)*rand+lambda/2];
Nmax=round((l-lambda)*(p-lambda)*(h-lambda)*.74/(4/3*pi*(3e8/f0/4)^3))
while (length(X)<1500)
    
    X_0=(l-lambda)*rand+lambda/2;
    Y_0=(p-lambda)*rand+lambda/2;
    Z_0=(h-lambda)*rand+lambda/2;
    
    D=sqrt((X_0-X).^2+(Y_0-Y).^2+(Z_0-Z).^2);
    if min(D)>lambda/2
        X=[X;X_0];
        Y=[Y;Y_0];
        Z=[Z;Z_0];
    end
end
for v=1:1:61
    disp([num2str(v),'/61']) 
tic
    filename = sprintf('%delem_%dnssourcesproches%d.mat',nbre_elements,round(Lt/(1e-9)),v);
    load(filename)
    for u=1:length(X)
        
        X_1=X(u);
        Y_1=Y(u);
        Z_1=Z(u);
        
       % disp('CIR')
        [Sx,Sy,Sz]=CIR(X_1,Y_1,Z_1);
        
        
        %Frequency response
        %disp('FFT')
        
        
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
        FFTx(v,u,:)=abs(Yx(1:4000));
        FFTy(v,u,:)=abs(Yy(1:4000));
        FFTz(v,u,:)=abs(Yz(1:4000));
        freq=f(1:4000);
    end
toc
end
save('Resultsourcesproches.mat','X','Y','Z','FFTx','FFTy','FFTz','freq','Nmax','f0')
