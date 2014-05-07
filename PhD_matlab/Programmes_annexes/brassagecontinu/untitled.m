clear all

load SR360_1000MHZ_360.mat
theta=1:360
SR=(Amplitude_crete/mean(Amplitude_crete)).^2;

thetai = .01:.01:360; 
SRi = interp1(theta,SR,thetai,'spline'); 





thetap=2*pi/60; %vitesse de rotation du brasseur en rad/s.
Lt=1/thetap*2*pi
T=2.3%periode de fonctionnement
rho_f=.5 %rapporty cyclique de fonctionnement de la machine
dt=1/thetap*2*pi/length(thetai);
t=dt:dt:Lt;
t1=dt:dt:T;%une période
Fu=ones(round(rho_f*length(t1)),1);
Fu=[Fu;zeros(length(t1)-length(Fu),1)];

Functionning=Fu;
while length(Functionning)<length(t)
    Functionning=[Functionning;Fu];
end
Functionning=Functionning(1:length(t));

power=SRi';
while length(power)<length(t)
    power=[power;SRi'];
end
power=power(1:length(t));

essai=power.*Functionning;
essaidB=10*log10(essai);

critere=0 %en dBm
U=find(essaidB>critere);

for u=1:length(U)
    Susceptibility(U(u))=essaidB(U(u));
end

Susceptibility=[Susceptibility zeros(1,length(t)-length(Susceptibility))];

figure(1)
hold on
plot(t,10*log10(power),'--b')
plot(t,essaidB,'-g','LineWidth',1)
plot(t,Susceptibility,'-r','LineWidth',2)
plot(t,critere*(ones(1,length(t))),'--k','LineWidth',2)

grid on
xlabel('Time [s]')
ylabel('Level [dB]')
hold off
pfailure=length(U)/length(essaidB)


%temps d'exposition max
for i=length(U):-1:2
    Z(i)=U(i)-U(i-1);
end

for j=1:length(Z)
    if Z(j)~=1
        Z(j)=0;
    end
end

ZZ=cumsum(Z);
h=0;
for i=2:length(ZZ)
    if ZZ(i-1)==ZZ(i)
        ttexp(i)=ZZ(i);
    else
        h=h+1;
    end
end
if h==length(ZZ)-1
    tmax=max(t)
    
else
    TEXP=find(ttexp>0);
    for i=length(TEXP):-1:2
        texp(i)=ttexp(TEXP(i))-ttexp(TEXP(i-1));
    end
    tmax=max(texp)*(t(2)-t(1))
end


