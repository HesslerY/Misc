alp=pi*rand(1,150);
bet=2*pi*rand(1,150);



D=1;% taille du syst�me en lambda
n =5;
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
frant = fres(thetaa, phia, n, xa, ya, za, Ia);
faant = fantres(frant, thetaa, phia, S);
d1ant = Directivite(faant, thetaa, phia)
tracediag(faant,thetaa, phia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



save('antenne12dB.mat','D','Ia','S','d1ant','dph','n','xa','ya','za')