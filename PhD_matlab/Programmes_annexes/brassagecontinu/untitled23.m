x=0:.1:23
y=log(.1)./log(1-exp(-x))

plot(x,y)


epsilon=exp(log(exp(-1*x)))
plot(x,epsilon*360*0.95)


lambda=1

%prob>c dB
k1=log(.5)/log(1-(exp(-lambda*10^(4/10))))

%prob = c dB

k2=log(0.1)/log(1-(exp(-lambda*10^(4/10))-(exp(-lambda*10^(6/10)))))

