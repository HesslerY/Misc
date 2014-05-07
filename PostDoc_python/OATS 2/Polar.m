clear java
clear all
c=2.998e8;


R = 10;  %radius (m)
h = 1;   %heigth of the EUT (m)
Reflec = 1; %Reflection coefficient
%signal
f=10e6:20e6:3e9;
lambda=c./f;

%EUT
%positions of the dipoles and amplitude and phase

%random EUT
L=.5; %(m)
n=20; %number of dipoles
for u=1:10000
x=L*rand(n,1)-L/2;
y=L*rand(n,1)-L/2;
z=L*rand(n,1)-L/2;
tilt=acos(2*rand(n,1)-1);
azimut=2*pi*rand(n,1);

amplitude=rand(n,1);
phas=2*pi*rand(n,1);

I = [x y z tilt azimut amplitude phas];

%%%OATS
I(:,3)=I(:,3)+h;
Ip = I;
Ip(:,3)=-Ip(:,3);
Ip(:,6)=Reflec*Ip(:,6);

POS=[I;Ip];

ph=2*pi*rand;
%Pac=[cos(ph)*I(:,6).*cos(I(:,7)).*cos(I(:,5)).*sin(I(:,4))+sin(ph)*I(:,6).*cos(I(:,7)).*sin(I(:,5)).*sin(I(:,4)) I(:,6).*cos(I(:,7)).*cos(I(:,4))];
%P=[cos(ph)*POS(:,6).*cos(POS(:,7)).*cos(POS(:,5)).*sin(POS(:,4)+sin(ph)*POS(:,6).*cos(POS(:,7)).*sin(POS(:,5)).*sin(POS(:,4)) POS(:,6).*cos(POS(:,7)).*cos(POS(:,4))];

Pac=[I(:,6).*cos(I(:,5)).*sin(I(:,4)) I(:,6).*cos(I(:,4))];
P=[POS(:,6).*cos(POS(:,5)).*sin(POS(:,4)) POS(:,6).*cos(POS(:,4))];



PPac(u,:)=abs(sum(Pac));
PP(u,:)=abs(sum(P));
end
mean(PPac)
mean(PP)