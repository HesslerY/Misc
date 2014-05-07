clear all

load CIR1.mat
load CIR2.mat
load CIR3.mat
load CIR4.mat

Sx=[Sx1;Sx2;Sx3;Sx4];
Sy=[Sy1;Sy2;Sy3;Sy4];
Sz=[Sz1;Sz2;Sz3;Sz4];

clear Sx1 Sx2 Sx3 Sx4 Sy1 Sy2 Sy3 Sy4 Sz1 Sz2 Sz3 Sz4


for i=1:150
h=hamming(N)';

Sxh=h.*Sx(i,:);
Syh=h.*Sy(i,:);
Szh=h.*Sz(i,:);

Sxf=Sx(i,:);
Syf=Sy(i,:);
Szf=Sz(i,:);
    
    
    
    
    
%FFT

Fs = N/Lt; %sampling frequency
T = 1/Fs;                     % sample length
L = N;                  %Number of points
tt = (0:L-1)*T;                % time...

NFFT = 2^nextpow2(N); % Next power of 2 from length of y
Yx = fft((Sxf),NFFT)/N;
Yy = fft((Syf),NFFT)/N;
Yz = fft((Szf),NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2);

% Three axis frequency response
FFTx(i,:) = abs(Yx(1:NFFT/2));
FFTy(i,:) = abs(Yy(1:NFFT/2));
FFTz(i,:) = abs(Yz(1:NFFT/2));
freq = f;

%FFT Hamming
Fs = N/Lt; %sampling frequency
T = 1/Fs;                     % sample length
L = N;                  %Number of points
tt = (0:L-1)*T;                % time...

NFFT = 2^nextpow2(N); % Next power of 2 from length of y
Yxh = fft((Sxh),NFFT)/N;
Yyh = fft((Syh),NFFT)/N;
Yzh = fft((Szh),NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2);

% Three axis frequency response
FFTxh(i,:) = abs(Yxh(1:NFFT/2));
FFTyh(i,:) = abs(Yyh(1:NFFT/2));
FFTzh(i,:) = abs(Yzh(1:NFFT/2));
freqh = f;

end


for i=1:16000
wblfit(FFTzh(:,1500))
    

end
    









