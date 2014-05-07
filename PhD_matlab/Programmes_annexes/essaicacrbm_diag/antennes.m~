alp=pi/6;
bet=pi/4;



D=1;% taille du syst�me en lambda
n =3;
dx=0.25;
xa = linspace(0,n*dx,n)
ya = ones(n,1);
za = zeros(n,1);
[Ia, dph] = dolph3(1, dx, n,10);
%Ia=ones(1,n)
dph=zeros(1,n-1);
S = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Diag
nt = 100;
np = 191;
theta2 = linspace(0, pi, nt);
phi2 = linspace(0, 2*pi, np);
[thetaa, phia] = meshgrid(theta2, phi2);
frant = fres(mod(thetaa-alp,pi), mod(phia-bet,2*pi), n, xa, ya, za, Ia);
faant = fantres(frant, mod(thetaa-alp,pi), mod(phia-bet,2*pi), S);
d1ant = Directivite(faant, thetaa, phia)
tracediag(faant,thetaa, phia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

10*log10(d1ant)

save('antenne7dB.mat','D','Ia','S','d1ant','dph','n','xa','ya','za')