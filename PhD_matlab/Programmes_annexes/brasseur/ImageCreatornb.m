%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                               %%%%%%
%%%%%                       IMAGE CREATOR N                         %%%%%%
%%%%%           by E. Amador (emmanuel.amador@gmail.com)            %%%%%%
%%%%%                         IETR/DGA                              %%%%%%
%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

Lt=1e-6; %Time-window length in seconds
f0=1e9;
c=3e8;
%Physical dimensions of the cavity (length, width, heigth) in meters
l=8.7;
p=3.7;
h=2.9;

filename=['1delem_1000nssourcesproches1x.mat';'1delem_1000nssourcesproches2x.mat';'1delem_1000nssourcesproches3x.mat';'1delem_1000nssourcesproches4x.mat';'1delem_1000nssourcesproches5x.mat';'1delem_1000nssourcesproches6x.mat';
'1delem_1000nssourcesproches1y.mat';'1delem_1000nssourcesproches2y.mat';'1delem_1000nssourcesproches3y.mat';'1delem_1000nssourcesproches4y.mat';'1delem_1000nssourcesproches5y.mat';'1delem_1000nssourcesproches6y.mat';
'1delem_1000nssourcesproches1z.mat';'1delem_1000nssourcesproches2z.mat';'1delem_1000nssourcesproches3z.mat';'1delem_1000nssourcesproches4z.mat';'1delem_1000nssourcesproches5z.mat';'1delem_1000nssourcesproches6z.mat'];


lambda=c/f0;
f0=1e9; %monochromatic pulse frequency

%Reception point coordinates

lambda=c/f0;

X=[2.2*rand+7.3];
Y=[2.2*rand+0.3];
Z=[(2.3)*rand+0.3];

Nmax=round(2.2*2.2*2.3*.74/(4/3*pi*(3e8/f0/4)^3))
while (length(X)<100)

    X_0=2.2*rand+7.3;
    Y_0=2.2*rand+0.3;
    Z_0=2.3*rand+0.3;

    D=sqrt((X_0-X).^2+(Y_0-Y).^2+(Z_0-Z).^2);
    if min(D)>lambda/2
        X=[X;X_0];
        Y=[Y;Y_0];
        Z=[Z;Z_0];
    end
end
X=[X;X;X];
Y=[Y;Y;Y];
Z=[Z;Z;Z];

amplitude=ones(300,1);
Phaseelement=zeros(300,1);
tilt=[pi/2*ones(200,1); zeros(100,1)];
azimut=[zeros(100,1); pi/2*ones(100,1);zeros(100,1)];

disp('matrice des source generee')
for iii=1:length(X)
tic
    disp([num2str(iii),'/',num2str(length(X))]) 
POS=[];
POSP=[];
POSI=[];


c=3e8;
dmax=Lt*c; %maximal distance for the choses time-window
ordre=round(dmax/min([l;p;h]))+1; %Maximum order

Memo=8*pi*dmax^3/l/p/h*80;
disp([num2str(round(Memo/1e6)),' MB needed'])

%Generation of the k=0 horizontal plane in four quadrants
disp('Even plane...')
for z=iii
    
    
    
    for i=1:ordre
        for j=1:ordre
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP=[POSP; 2*i*l-X(z)	2*j*p+Y(z)	Z(z)    abs(2*i)-1  abs(2*j)      0	Phaseelement(z)	amplitude(z)   pi-tilt(z) mod(2*pi-azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p+Y(z)	Z(z)	abs(2*i)    abs(2*j)      0	Phaseelement(z)	amplitude(z)   tilt(z) mod(azimut(z),2*pi);
                    2*i*l-X(z)	2*j*p-Y(z)	Z(z)    abs(2*i)-1  abs(2*j)-1    0	Phaseelement(z)	amplitude(z)   tilt(z)   mod(pi+azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p-Y(z)	Z(z)	abs(2*i)    abs(2*j)-1	  0	Phaseelement(z)	amplitude(z)   pi-tilt(z)  mod(pi-azimut(z),2*pi)];
            end
        end
    end
    
    
    
    for i=-ordre:0
        for j=1:ordre
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP=[POSP; 2*i*l-X(z)	2*j*p+Y(z)	Z(z)    abs(2*i)+1  abs(2*j)      0	Phaseelement(z)	amplitude(z)   pi-tilt(z) mod(2*pi-azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p+Y(z)	Z(z)	abs(2*i)    abs(2*j)      0	Phaseelement(z)	amplitude(z)    tilt(z) mod(azimut(z),2*pi);
                    2*i*l-X(z)	2*j*p-Y(z)	Z(z)    abs(2*i)+1  abs(2*j)-1    0	Phaseelement(z)	amplitude(z)    tilt(z)  mod(pi+azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p-Y(z)	Z(z)	abs(2*i)    abs(2*j)-1	  0	Phaseelement(z)	amplitude(z)   pi-tilt(z)  mod(pi-azimut(z),2*pi)];
            end
        end
    end
    
    
    
    for i=-ordre:0
        for j=-ordre:0
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP=[POSP; 2*i*l-X(z)	2*j*p+Y(z)	Z(z)    abs(2*i)+1  abs(2*j)      0	Phaseelement(z)	amplitude(z)   pi-tilt(z) mod(2*pi-azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p+Y(z)	Z(z)	abs(2*i)    abs(2*j)      0	Phaseelement(z)	amplitude(z)    tilt(z) mod(azimut(z),2*pi);
                    2*i*l-X(z)	2*j*p-Y(z)	Z(z)    abs(2*i)+1  abs(2*j)+1    0	Phaseelement(z)	amplitude(z)    tilt(z)  mod(pi+azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p-Y(z)	Z(z)	abs(2*i)    abs(2*j)+1	  0	Phaseelement(z)	amplitude(z)   pi-tilt(z)  mod(pi-azimut(z),2*pi)];
            end
        end
    end
    
    
    for i=1:ordre
        for j=-ordre:0
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP=[POSP; 2*i*l-X(z)	2*j*p+Y(z)	Z(z)    abs(2*i)-1  abs(2*j)      0	Phaseelement(z)	amplitude(z)   pi-tilt(z) mod(2*pi-azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p+Y(z)	Z(z)	abs(2*i)    abs(2*j)      0	Phaseelement(z)	amplitude(z)    tilt(z) mod(azimut(z),2*pi);
                    2*i*l-X(z)	2*j*p-Y(z)	Z(z)    abs(2*i)-1  abs(2*j)+1    0	Phaseelement(z)	amplitude(z)    tilt(z)  mod(pi+azimut(z),2*pi);
                    2*i*l+X(z)	2*j*p-Y(z)	Z(z)	abs(2*i)    abs(2*j)+1	  0	Phaseelement(z)	amplitude(z)   pi-tilt(z)  mod(pi-azimut(z),2*pi)];
            end
        end
    end
    
    
    
    
    
    toc
end

%duplication along the vertical axis

disp('Odd plane...')
%odd horizontal plane generation 
POSI=POSP; 
POSI(:,3)=h-POSI(:,3);

POSI(:,10)=mod(POSI(:,10)+pi,2*pi);%inversion de l'azimuth

disp('Oz axis duplication')
for k=-ordre:ordre
    if mod(k,2)==0
        POS=[POS;POSP+[zeros(length(POSP),1) zeros(length(POSP),1) (k*h)*ones(length(POSP),1)   zeros(length(POSP),1)   zeros(length(POSP),1)  abs(k)*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
    else
        POS=[POS;POSI+[zeros(length(POSI),1) zeros(length(POSI),1) (k*h)*ones(length(POSI),1)   zeros(length(POSP),1)   zeros(length(POSI),1)  abs(k)*ones(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1)]];
    end
    
end
toc
disp('Saving...')
%saving the POSITION matrix
%filename = sprintf('1delem_%dnssourcesproches%d.mat',round(Lt/(1e-9)),iii);
if iii<101
    filename=sprintf('x%d.mat',iii)
else
    if iii<201
        filename=sprintf('y%d.mat',mod(iii,100)+1)

    else
        filename=sprintf('z%d.mat',mod(iii,100)+1)
    end
end
save(filename,'POS')
toc
disp('Done.')
end
save('sources.mat','X','Y','Z')


