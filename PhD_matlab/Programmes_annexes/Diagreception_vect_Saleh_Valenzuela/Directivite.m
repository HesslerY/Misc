

function D = Directivite(Fa, theta, phi)

nt = length(theta(1,:));
np = length(phi(:,1));

dtheta = pi/nt;
dphi = (2*pi)/np;

%Calcul de l'angle solide
omega = sum(sum((Fa .^ 2) .* sin(theta)*dtheta*dphi));

D = 4*pi / omega;
