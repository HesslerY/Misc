clear

p=10000;
t=(1:p)/p;
s=sin((4*t).^2*pi);
tt=t.^2;
ss=s;

q=1;
t2=1/(p*q):1/(p*q):max(sqrt(tt));
tt2=(2/(p*q)):(2/(p*q)):max((tt));
s2=zeros(1,length(t2));
s3=zeros(1,length(tt2));
 for i=1:length(tt)
     if (round(length(t2)*sqrt(tt(i))/max(sqrt(tt)))+1)<length(t2)
     s2(round(length(t2)*sqrt(tt(i))/max(sqrt(tt)))+1)=ss(i);
         
     end
 end
 
 for i=1:length(tt)
     if (round(length(tt2)*(tt(i))/max((tt)))+1)<length(tt2)
     s3(round(length(tt2)*(tt(i))/max((tt)))+1)=ss(i);
         
     end
 end


%Frequency response
disp('FFT')


Fs = length(tt)/max(t2); %sampling frequency
T = 1/Fs;                     % sample length
LL = length(t2);                  %Number of points

NFFT = 2^nextpow2(LL); % Next power of 2 from length of y
Yz = fft((s2),NFFT);
f = Fs/2*linspace(0,1,NFFT/2);
FFTz=abs(Yz(1:NFFT/2));
freq=f;


 

 
 %Frequency response
disp('FFT')


Fst = length(tt)/max(tt); %sampling frequency
Tt = 1/Fst;                     % sample length
Lt = length(tt);                  %Number of points

NFFTt = 2^nextpow2(Lt); % Next power of 2 from length of y
Yzt = fft((s3),NFFTt);
ft = Fst/2*linspace(0,1,NFFTt/2);
FFTzt=abs(Yzt(1:NFFTt/2));
freqt=ft;



figure(2)
plot(freq,10*log10(FFTz),(freqt),10*log10(FFTzt))
