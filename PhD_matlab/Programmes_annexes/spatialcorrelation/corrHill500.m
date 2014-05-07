%reproduction resultat de Hill


clear

global Lt c Rx Ry Rz N POS

freq=[];
l=8.7;
p=3.7;
h=2.9;

c=3e8;%
Lt=3e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance

tau=.3e-6; %length of the pulse in seconds
f0=1e9; %monochromatic pulse frequency

N=3e6%round(30*Lt*f0) %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.998;
Ry=0.998;
Rz=0.998;


%Reception point coordinates

lambda=c/f0;

X=[.8]
Y=[.8]
Z=[.8]

M=500
while(length(X)<M)
    a=length(X);
    pr=(8.7+2.9+3.7)*rand;
    if pr<2.9
        X=[X;X(a)];
        Y=[Y;Y(a)+0.015];
        Z=[Z;Z(a)];
    else
        if pr<(2.9+2)
            X=[X;X(a)];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)+0.015];
        else

            X=[X;X(a)+0.015];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)];
        end
    end
end
max([X Y Z])

X=[X;1.6];
Y=[Y;1];
Z=[Z;1];

while(length(X)<2*M)
    a=length(X);
    pr=(8.7+2.9+3.7)*rand;
    if pr<2.9
        X=[X;X(a)];
        Y=[Y;Y(a)+0.015];
        Z=[Z;Z(a)];
    else
        if pr<(2.9+2)
            X=[X;X(a)];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)+0.015];
        else

            X=[X;X(a)+0.015];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)];
        end
    end
end
max([X Y Z])

X=[X;2.4];
Y=[Y;1];
Z=[Z;1];

while(length(X)<3*M)
    a=length(X);
    pr=(8.7+2.9+3.7)*rand;
    if pr<2.9
        X=[X;X(a)];
        Y=[Y;Y(a)+0.015];
        Z=[Z;Z(a)];
    else
        if pr<(2.9+2)
            X=[X;X(a)];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)+0.015];
        else

            X=[X;X(a)+0.015];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)];
        end
    end
end
max([X Y Z])



load('1elemIETR_3000ns.mat')
for u=1:length(X)
tic
disp([num2str(u),'/1500'])
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
    FFTx(u,:)=abs(Yx(1:200000));
    FFTy(u,:)=abs(Yy(1:200000));
    FFTz(u,:)=abs(Yz(1:200000));
    freq=f(1:200000);
    toc
end
toc
save('ResultHill3x500elements.mat','X','Y','Z','FFTx','FFTy','FFTz','freq','f0')
