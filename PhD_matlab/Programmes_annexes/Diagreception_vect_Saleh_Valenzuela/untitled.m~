n =2;
dx=5
xa = 1:dx:dx*n+1;%D*rand(n,1);
ya = ones(1,n);%D*rand(n,1);
za = ones(1,n);%D*rand(n,1);

%[Ia, dph] = dolph3(1, dx, n,10)
Ia=ones(1,n)
dph=zeros(1,n-1)
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
d1ant = Directivite(faant, thetaa, phia);
tracediag(faant,thetaa, phia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=10*log10(d1ant)