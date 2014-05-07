function bide = tracediag(Fa,varargin)

% Fa(nphi+1,ntheta+1) est la matrice de la fonction caractéristique normalisée
% évaluée aux divers angles d'observation theta et phi
% [Theta Phi] = meshgrid(theta,phi)
%
% utilise la fonction polair definie par DG (dans outil)
% qui permet de choisir le rayon maximum via un paramètre
% le sens de rotation et l'angle de l'axe x
%
% tracediag(Fa)
% tracediag(Fa,Theta,Phi)
% tracediag(Fa,'PARAM',valeur,...)
% tracediag(Fa,Theta,Phi,'PARAM',valeur...)
% 'lin' ou 'log' permet de choisir l'échelle linéaire ou dB
% 'max',valeur fixe l'écart maximal des données affichées (axe radial)
% 'azimut',angle et 'elevation',angle donnent l'angle de vision pour le graphe 3D
% 'theta',[angle,angle] 'phi',[angle,angle] limitent la portion du graphe 3D
% 'xy', 'xz', 'yz' et/ou 'xyz' sélectionnent le ou les graphes désirés

d2r = pi/180;
line = 'b-';
 
ntheta = size(Fa,2)-1;
nphi   = size(Fa,1)-1;
if (nargin==1) | isstr(varargin{1})
	dtheta = pi/ntheta;
	dphi   = 2*pi/nphi;
	theta  = 0:dtheta:pi;
	thetaa = 0:dtheta:2*pi;
    phi    = 0:dphi:2*pi;
    [Theta Phi] = meshgrid(theta,phi);
    
	Faxy = Fa(:,ntheta/2+1)';
	Faxz = [Fa(1,:) Fa(fix(nphi/2)+1,end-1:-1:1)];
    Fayz = [Fa(fix(nphi/4)+1,:) Fa(fix(3*nphi/4)+1,end-1:-1:1)];
    
    nin = 1;	
else
    Theta  = varargin{1};
    Phi    = varargin{2};
    dtheta = Theta(1,2)-Theta(1,1);
    dphi   = Phi(2,1)-Phi(1,1);
    theta  = Theta(1,:);
    thetaa = [theta 2*pi-theta(end-1:-1:1)];
    phi    = Phi(:,1)';
    
    [tmp ind]  = min(abs(Theta'-pi/2));
    Faxy = (tmp(1)<dtheta)*Fa(:,ind(1))';
	[tmpa inda]  = min(abs(Phi));
    [tmpb indb]  = min(abs(Phi-pi));
    Faxz = [(tmpa(1)<dphi)*Fa(inda(1),:) (tmpb(1)<dphi)*Fa(indb(1),end-1:-1:1)];
	[tmpa inda]  = min(abs(Phi-pi/2));
	[tmpb indb]  = min(abs(Phi-3*pi/2));
    Fayz = [(tmpa(1)<dphi)*Fa(inda(1),:) (tmpb(1)<dphi)*Fa(indb(1),end-1:-1:1)];
    
    nin = 3;
end

%valeurs initiales
lindBn = 1;
maxe = 0;
azimut = -36;
elevat = 30;
grafs = [0 0 0 0];
thlim = [0 pi];
phlim = [0 2*pi];

while nin<=nargin-1
	switch varargin{nin}
    case 'lin'
        lindBn = 1;
        nin = nin+1;
    case 'log'
        lindBn = 0;
        nin = nin+1;
    case 'max'
        maxe = varargin{nin+1};
        nin = nin+2;
    case 'azimut'
        azimut = varargin{nin+1};
        nin = nin+2;
    case 'elevation'
        elevat = varargin{nin+1};
        nin = nin+2;
    case 'theta'
        thlim = varargin{nin+1}*d2r;
        if thlim(1)<0 error('borne theta négative'); end
        nin = nin+2;
     case 'phi'
        phlim = varargin{nin+1}*d2r;
        if phlim(1)<0 phlim(1) = phlim(1)+2*pi; end
        nin = nin+2;
   case 'xy'
        grafs(1) = 1;
        nin =nin+1;
    case 'xz'
        grafs(2) = 1;
        nin = nin+1;
    case 'yz'
        grafs(3) = 1;
        nin = nin+1;
    case 'xyz'
        grafs(4) = 1;
        nin = nin+1;
    otherwise
        %error('mauvais paramètre');
        line = varargin{nin};
        nin = nin+1;
    end
end

if (sum(grafs)==0) grafs = [1 1 1 1]; end
if lindBn
    if ~maxe maxe = 1; end
	text = 'norm';
else %log
    if ~maxe maxe = 40; end
    Fa   = maxe+20*log10(max(Fa,10^(-maxe/20)));
    Faxy = maxe+20*log10(max(Faxy,10^(-maxe/20)));
    Faxz = maxe+20*log10(max(Faxz,10^(-maxe/20)));
    Fayz = maxe+20*log10(max(Fayz,10^(-maxe/20)));
    text = 'dB';
end;

ngpl = sum(grafs(1:3));
if ngpl
    figure;
    i = 1;
    if grafs(1)
        subplot(1,ngpl,i), polair(phi,Faxy,maxe,90,1,line);
        xlabel(['Plan xy ' text ' (\phi,\theta=90)']), axis(maxe*[-1.15 1.15 -1 1]);
        i = i+1;
    end
    if grafs(2)
        subplot(1,ngpl,i), polair(thetaa,Faxz,maxe,270,1,line);
        xlabel(['Plan xz ' text ' (\phi=0,\theta)']),  axis(maxe*[-1.15 1.15 -1 1]);
        i = i+1;
    end
    if grafs(3)

        subplot(1,ngpl,i), polair(thetaa,Fayz,maxe,90,-1,line);
    	xlabel(['Plan yz ' text ' (\phi=90,\theta)']), axis(maxe*[-1.15 1.15 -1 1]);
    end
end

if grafs(4)
	it = find((theta>=thlim(1)) & (theta<=thlim(2)));
	if phlim(1)<phlim(2)
    	ip = find((phi>=phlim(1)) & (phi<=phlim(2)));
	else
    	ip = find((phi>=phlim(1)) | (phi<=phlim(2)));
    end
    Fa = Fa(ip,it);
    Theta = Theta(ip,it); Phi = Phi(ip,it);
    
	figure;
	view(azimut,elevat);
	hidden off;
	su=mesh(Fa.*cos(Phi).*sin(Theta),Fa.*sin(Phi).*sin(Theta),Fa.*cos(Theta),Fa);
	title(['3D ' text]), axis(maxe*[-1 1 -1 1 -1 1]);
	cmap=hot;
	set(su,'FaceColor','interp','EdgeColor','none');
%	set(su,'FaceColor','interp');
	set(su,'AmbientStrength',.6);
%	colormap(cmap(6:end-3,[3 2 1])); % de bleu foncé à bleu pâle
    colormap(cmap(end-3:-1:6,[1 2 3])); % de jaune à rouge foncé
end
