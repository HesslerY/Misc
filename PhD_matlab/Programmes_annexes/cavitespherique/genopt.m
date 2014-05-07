R=.15;
X0=0;
Y0=0;
Z0=0;

amplitude=1
phaselement=1
tilt=0
azimut=0
theta=0:2*pi/200:pi;
phi=0:2*pi/200:2*pi;
POS=[0 0 0 0 1 0 0 0]
POSI=[];
POSP=[];
ordre=200;

for j=1:1:length(theta)
    for k=1:1:length(phi)
        POSP=[POSP;2*R*cos(phi(k))*sin(theta(j)) 2*R*sin(phi(k))*sin(theta(j)) 2*R*cos(theta(j)) 0 1 0 0 0];
        POSI=[POSI;2*R*cos(phi(k))*sin(theta(j)) 2*R*sin(phi(k))*sin(theta(j)) 2*R*cos(theta(j)) 0 1 0 rem(2*theta(j),2*pi) phi(k)];
    end
end
for i=1:1:ordre
    if mod(i,2)==0
        
                POS=[POS;POSP+[(i-1)*POSP(:,1) (i-1)*POSP(:,2) (i-1)*POSP(:,3) i*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
        
    else
        
                POS=[POS;POSI+[(i-1)*POSI(:,1) (i-1)*POSI(:,2) (i-1)*POSI(:,3) i*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
       
    end
end
save sphere200.mat 'POS'
FFTplan


