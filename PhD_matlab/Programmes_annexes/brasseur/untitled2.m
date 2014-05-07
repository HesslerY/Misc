clear

p=1000;
t=(1:p)/p;
s=sin((4*t).^2*pi);

M=5;

tt=t.^2;
ss=s;
for i=1:M
    tt=[tt t.^2+max(tt)];
    ss=[ss s];
end


q=1;
t2=1/(p*q):1/(p*q):max(sqrt(tt));
s2=zeros(1,length(t2));
 for i=1:length(tt)
     if (round(length(t2)*sqrt(tt(i))/max(sqrt(tt)))+1)<length(t2)
     s2(round(length(t2)*sqrt(tt(i))/max(sqrt(tt)))+1)=ss(i);
         
     end
 end
%s2=interpft(ss,length(t2))

figure(1)
subplot(3,2,1)
plot(t,s,'.','MarkerSize',2)
xlim([0 max(t)])
xlabel('$t$','Interpreter','Latex')
title('Original signal vs. $t$','Interpreter','Latex')
subplot(3,2,2)
plot(t.^2,s,'.','MarkerSize',2)
xlim([0 max(t.^2)])
title('Original signal vs. $t^2$','Interpreter','Latex')
xlabel('$t^2$','Interpreter','Latex')
subplot(3,2,3)
plot(tt,ss,'.','MarkerSize',2)
xlim([0 max(tt)])
title('Extended signal vs. $t^2$','Interpreter','Latex')
xlabel('$t^2$','Interpreter','Latex')
subplot(3,2,4)
plot(sqrt(tt),ss,'.','MarkerSize',2)
xlim([0 max(sqrt(tt))])
title('Extended signal vs. $t$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')

subplot(3,2,[5 6])
plot(t2,s2,'.','MarkerSize',2)
xlim([0 max(t2)])
title('Extended signal vs. $t$ resampled','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')


%Frequency response
disp('FFT')


Fs = length(t)/max(t); %sampling frequency
T = 1/Fs;                     % sample length
L = length(t);                  %Number of points

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Yz = fft((s),NFFT);
f = Fs/2*linspace(0,1,NFFT/2);
FFTz=abs(Yz(1:NFFT/2));
freq=f;


%Frequency response
disp('FFT')


Fs = length(tt)/max(t2); %sampling frequency
T = 1/Fs;                     % sample length
LL = length(t2);                  %Number of points

NFFT = 2^nextpow2(LL); % Next power of 2 from length of y
Yz = fft((s2),NFFT);
f = Fs/2*linspace(0,1,NFFT/2);
FFTzp=abs(Yz(1:NFFT/2));
freqp=f;

figure(2)
plot(freq,10*log10(FFTz),freqp,10*log10(FFTzp))
