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

N=50;
   
Lt=1e-6; %Time-window length in seconds
f0=.5e9;
c=3e8;
lambda=c/f0;

%Physical dimensions of the cavity (length, width, heigth) in meters
l=8.7;
p=3.7;
h=2.9;


XX=[(l-lambda)*rand+lambda/2];
YY=[(p-lambda)*rand+lambda/2];
ZZ=[(h-lambda)*rand+lambda/2];

while (length(XX)<N)
    
    X_1=(l-lambda)*rand+lambda/2;
    Y_1=(p-lambda)*rand+lambda/2;
    Z_1=(h-lambda)*rand+lambda/2;
    
    D=sqrt((X_1-XX).^2+(Y_1-YY).^2+(Z_1-ZZ).^2);
    if min(D)>lambda/2
        XX=[XX;X_1];
        YY=[YY;Y_1];
        ZZ=[ZZ;Z_1];
       
    end

end
disp('matrice des source generee')
for iii=1:N
tic
    disp([num2str(iii),'/50']) 
    X=XX(iii);
    Y=YY(iii);
    Z=ZZ(iii);
%Radiating element(s)
%position

%amplitude
amplitude=[1];
%phase for array antennas synthesis or MIMO systems
Phaseelement=[0];

%angular orientation
tilt=[pi/2-acos(sqrt(2/3))];
azimut=[pi/4];

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
for z=1:1:length(X)
    
    
    
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
filename = sprintf('%delem_%dns%d.mat',length(X),round(Lt/(1e-9)),iii);
save(filename,'POS')
toc
disp('Done.')
end



