R=3;
X0=0;
Y0=0;
Z0=0;

amplitude=1
phaselement=1
tilt=pi/2
azimut=0
theta=0:2*pi/100:2*pi;
POS=[0 0 0 0 1 0 pi/2 0]
ordre=200;
for i=1:1:ordre
    if mod(i,2)==0
        for j=1:1:length(theta)
            POS=[POS;2*R*i*cos(theta(j)) 2*R*i*sin(theta(j)) 0 i 1 0 pi/2 0];
        end
    else
        for j=1:1:length(theta)
            POS=[POS;2*R*i*cos(theta(j)) 2*R*i*sin(theta(j)) 0 i 1 0 pi/2 rem(2*theta(j),2*pi)];
        end
    end
end



