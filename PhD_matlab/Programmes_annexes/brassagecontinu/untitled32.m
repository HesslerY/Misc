load SR360times25_300MHZ_LogP_9000.mat

SR300=Amplitude_crete;
SR300=SR300.^2/mean(SR300.^2);
load SR360times25_500MHZ_LogP_9000.mat

SR500=Amplitude_crete;
SR500=SR500.^2/mean(SR500.^2);
load SR360times25_1000MHZ_cornetcornet_9000.mat

SR1000=Amplitude_crete;
SR1000=SR1000.^2/mean(SR1000.^2);

theta=360/9000:360/9000:360

plot(theta,10*log10(SR1000))%,theta,SR500,theta,SR1000)
ylim([4 6])
load SRopldiff.mat
SR=SR';


i=410;
SR1=SR(:,i).^2/mean(SR(:,i).^2);
thetai=360/90000:360/90000:360;
plot(thetai,10*log10(SR1))
ylim([4 6])