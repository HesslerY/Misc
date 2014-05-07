clear all
close all
load SR360_1000MHZ_360.mat
SR1=(Amplitude_crete/mean(Amplitude_crete)).^2;
load SR360_3000MHZ_360.mat
theta=1:360;
SR2=(Amplitude_crete/mean(Amplitude_crete)).^2;
S0=sin(theta/16.7).^4;
SR0=(S0/mean(S0)).^2;



powerdB0=10*log10(SR0);
powerdB1=10*log10(SR1);
powerdB2=10*log10(SR2);

subplot(3,1,1)
plot(theta,powerdB0)
xlabel('$\theta$')
ylabel('dBm')
grid on
title('sin')
xlim([0 360])
ylim([0 max(powerdB0)])

subplot(3,1,2)
plot(theta,powerdB1)
xlabel('$\theta$')
ylabel('dBm')
grid on
title('1 GHz')
xlim([0 360])
ylim([0 max(powerdB1)])


subplot(3,1,3)
plot(theta,powerdB2)
xlabel('$\theta$')
ylabel('dBm')
grid on
xlim([0 360])
ylim([0 max(powerdB2)])
title('3 GHz')