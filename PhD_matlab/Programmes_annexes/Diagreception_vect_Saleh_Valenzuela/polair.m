function hpol = polair(theta,rho,maxrho,th0,sens,line)
% original /gel/logiciels/matlab5.1/toolbox/graph2d
% maxrho = rayon maximal fixé
% th0 = l'angle correpondant à l'axe x en degrés
% sens=1 : antihoraire (par défaut); sens=0 : horaire; sens=2 : 0-180 degrés dans les 2 sens

line_style = 'auto';

if nargin < 2
    error('Il faut 2 arguments d''entrée ou plus.')
elseif nargin == 2
    maxrho = max(abs(rho(:)));
    th0 = 0;
    sens = 1;
elseif nargin == 3
    th0 = 0;
    sens = 1;
elseif nargin == 4
    sens = 1;
end
ss = 1;
if sens == 0
    sens = -1;
elseif sens == 2
    sens = 1; ss= -1;
end
if nargin ~= 6
    line = 'b-';
end

if ~isequal(size(theta),size(rho))
    error('THETA et RHO doivent avoir la même dimension.');
end

[mr,nr] = size(rho);

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    axis image; v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
ticks = 5;
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'));
    end
% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1);
%        text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
        text((i-rinc/2),0, ...
            ['  ' num2str(i)],'verticalalignment','bottom','FontSize',7)
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12 - pi*th0/180;
    cst = cos(sens*th); snt = sin(sens*th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',1)

% annotate spokes in degrees
    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),'horizontalalignment','center','FontSize',7)
        if i == length(th)
            loc = int2str(0);
        else
           loc = int2str(180+ss*i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center','FontSize',7)
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
theta = sens*(theta - pi*th0/180);
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
%if strcmp(line_style,'auto')
%    q = plot(xx,yy);
%else
    q = plot(xx,yy,line);
%end
if nargout > 0
    hpol = q;
end
if ~hold_state
    axis image; axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')
