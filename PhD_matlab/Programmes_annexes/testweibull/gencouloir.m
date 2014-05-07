
tic
Lt=3e-6 %temps de la fen???tre temporelle

%TH=0:2*pi/35:2*pi
% for ttt=1:length(TH)
% 
% 
%     X=[];
%     Y=[];
%     Z=[];
%     Phaseelement=[];
%     amplitude=[];
% 
%     m=3;
%     %Matrice des noeuds
%     X0=[7.3 7.3 7.3];
%     Y0=[1.4 1.4 1.4];
%     Z0=[.2 1 1.8];
%     T=[0 2.2*pi/3 4*pi/3];
%     P=[pi/6 pi/3 pi/4];
%     n=3;
%     dphi=0.9*pi+pi/n;
%     %dphi=0
%     f0=1e9
%     d=0.2*(1-1/n)
%
%     j=1
%     for i=1:1:n
%         X=[X;X0(j)-n*d/2+d*i*cos(TH(ttt)+T(j))*cos(P(j))];
%         Y=[Y;Y0(j)-n/2*d+d*i*sin(TH(ttt)+T(j))*cos(P(j))];
%         Z=[Z;Z0(j)-n/2*d+d*i*sin(P(j))];
%         amplitude=[amplitude;1/(n*m)];
%         Phaseelement=[Phaseelement;i*dphi];
%     end

%     j=2
%     for i=1:1:n
%         X=[X;X0(j)-n*d/2+d*i*cos(TH(ttt)+T(j))*cos(P(j))];
%         Y=[Y;Y0(j)-n/2*d+d*i*sin(TH(ttt)+T(j))*cos(P(j))];
%         Z=[Z;Z0(j)-n/2*d+d*i*sin(P(j))];
%         amplitude=[amplitude;1/(n*m)];
%         Phaseelement=[Phaseelement;i*dphi];
%     end

%     j=3
%     for i=1:1:n
%         X=[X;X0(j)-n*d/2+d*i*cos(TH(ttt)+T(j))*cos(P(j))];
%         Y=[Y;Y0(j)-n/2*d+d*i*sin(TH(ttt)+T(j))*cos(P(j))];
%         Z=[Z;Z0(j)-n/2*d+d*i*sin(P(j))];
%         amplitude=[amplitude;1/(n*m)];
%         Phaseelement=[Phaseelement;i*dphi];
%     end
X=[1];
Y=[1];
Z=[2];
amplitude=[1];
Phaseelement=[0];
tilt=[pi/4];
azimut=[pi/4];

POS=[];
POSP=[];
POSI=[];
%dimensions de la cavit???
l=23
p=2;
h=2;

%calcul de l'ordre correspondant
c=3e8;
dmax=Lt*c %distance maximale des sources
ordre=round(dmax/min([l;p;h]))+1 %ordre maximal corrspondant
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
    POSI=POSP;
    POSI(:,3)=h-POSI(:,3);
    %POSI(:,9)=pi+POSI(:,9);%changement de' tilt
    POSI(:,10)=mod(POSI(:,10)+pi,2*pi);%inversion de l'azimuth
for k=-ordre:ordre
    if mod(k,2)==0
        POS=[POS;POSP+[zeros(length(POSP),1) zeros(length(POSP),1) (k*h)*ones(length(POSP),1)   zeros(length(POSP),1)   zeros(length(POSP),1)  abs(k)*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
    else
        POS=[POS;POSI+[zeros(length(POSI),1) zeros(length(POSI),1) (k*h)*ones(length(POSI),1)   zeros(length(POSP),1)   zeros(length(POSI),1)  abs(k)*ones(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1)]];
    end
    
end 
    
filename = sprintf('%delem_%dnscouloir111.mat',length(X),round(Lt/(1e-9)));   
save(filename,'POS')
toc
%clear POS

% 
% %%special antennes
% POSANT=[];
% for g=1:1:length(POS)
%     if (POS(g,4)+POS(g,5)+POS(g,6))<3 
%         ;
%         POSANT=[POSANT;POS(g,:)];
%         %disp(g)
%     end
% end
% 
% disp(POSANT)



