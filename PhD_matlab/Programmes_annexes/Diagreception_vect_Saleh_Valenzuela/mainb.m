%diag al�atoire

%function main()

clear all;
close all;

nt = 100;
np = 191;

theta2 = linspace(0, pi, nt);
phi2 = linspace(0, 2*pi, np);
[theta, phi] = meshgrid(theta2, phi2);

D=10% taille du syst�me en lambda

%TELECHARGER LE PROG DOLPH

n =100;
dx=0.25
x = D*rand(n,1);%1;%1:dx:dx*n+1;%
y = D*rand(n,1);%1;%ones(1,n);%
z = D*rand(n,1);%1;%ones(1,n);%
%I=1;
%[I, dph] = dolph3(1, dx, n, 15)
I=complex(rand(n,1),rand(n,1));
%I = [0.9628+1i*0.2701;-4.1837-1i*0.8703;9.8866+1i*1.3603;-15.8051-1i*1.0822;18.3655;-15.8051+1i*1.0822;9.8866-1i*1.3603;-4.1837+1i*0.8703;0.9628-1i*0.2701]

S = 2;

fr = fres(theta, phi, n, x, y, z, I);

fa = fantres(fr, theta, phi, S);

d1 = Directivite(fa, theta, phi);

10*log10(d1)


tracediag(fa,theta, phi);
