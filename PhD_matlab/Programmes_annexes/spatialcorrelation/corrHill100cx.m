%reproduction resultat de Hill


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
f0=1e9; %monochromatic pulse frequency

N=1e6%round(30*Lt*f0) %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.998;
Ry=0.998;
Rz=0.998;


%Reception point coordinates

lambda=c/f0;


X_x=[(l-lambda)*rand+lambda/2];
Y_y=[(p-lambda)*rand+lambda/2];
Z_z=[(h-lambda)*rand+lambda/2];

X=[X_x];
Y=[Y_y];
Z=[Z_z];

        for i=1:59
            X=[X;X_x+i*0.005];
            Y=[Y;Y_y];
            Z=[Z;Z_z];
        end
Nmax=round((l-lambda)*(p-lambda)*(h-lambda)*.74/(4/3*pi*(3e8/f0/4)^3))
while (length(X)<6000)

    X_0=(l-lambda)*rand+lambda/2;
    Y_0=(p-lambda)*rand+lambda/2;
    Z_0=(h-lambda)*rand+lambda/2;

    D=sqrt((X_0-X).^2+(Y_0-Y).^2+(Z_0-Z).^2);
    if min(D)>lambda/2

        X=[X;X_0];
        Y=[Y;Y_0];
        Z=[Z;Z_0];
        for i=1:60
            X=[X;X_0+i*0.005];
            Y=[Y;Y_0];
            Z=[Z;Z_0];
        end
    end
end

load('1elem_1000ns.mat')
for u=1:length(X)
tic
disp([num2str(u),'/',num2str(length(X))])
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
   xFFTx(u,:)=abs(Yx(1:10000));
   xFFTy(u,:)=abs(Yy(1:10000));
   xFFTz(u,:)=abs(Yz(1:10000));
   freq=f(1:10000);
   toc
end

save('ResultHill500elements_x.mat','X','Y','Z','xFFTx','xFFTy','xFFTz','freq','f0')
toc