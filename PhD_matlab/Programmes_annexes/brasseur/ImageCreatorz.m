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
%clear

Lt=1e-6; %Time-window length in seconds
f0=1e9;
c=3e8;
%Physical dimensions of the cavity (length, width, heigth) in meters
l=8.7;
p=3.7;
h=2.9;

r=1.2;
 %u=1%:6
    Z=[];
    X=[];
    Y=[];
    
    %theta0=[0; pi/3; 3*pi/4;1.1*pi;4*pi/3;7*pi/4]
    for v=1:360
        X0=7.2+r*cos(v*2*pi/360)
        Y0=1.4+r*sin(v*2*pi/360)
        Z0=(u-1)*.4+.5
        
        X=[X0];
        Y=[Y0];
        Z=[Z0];
        
        amplitude=[1];
        Phaseelement=[0];
        tilt=[0];
        azimut=[0];
        
        disp('matrice des source generee')
        for iii=1
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
            end
            filename=sprintf('z_%d_%d.mat',u,v)
            
            save(filename,'POS')
            toc
            disp('Done.')
        
        %save('sources.mat','X','Y','Z')
    end
%end
