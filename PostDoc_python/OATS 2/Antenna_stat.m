clear java
clear all
c=2.998e8;

M=50;
MC=500;
%open area test setup
Hmax = 4;%maximum heigth of the antenna (m)
Hmin = 1;%minimum heigth of the antenna (m)
R = 10;  %radius (m)
h = 1;   %heigth of the EUT (m)
Reflec = 1; %Reflection coefficient
%signal
f=1e6:5e6:2e9;
lambda=c./f;

Ephac=zeros(M,MC,length(f));
Ephacgp=zeros(M,MC,length(f));
Ephoats=zeros(M,MC,length(f));
 
Erac=zeros(M,MC,length(f));
Eracgp=zeros(M,MC,length(f));
Eroats=zeros(M,MC,length(f));
Ethac=zeros(M,MC,length(f));
Ethacgp=zeros(M,MC,length(f));
Ezoats=zeros(M,MC,length(f));


%EUT
%positions of the dipoles and amplitude and phase
tic
for o=1:MC
    %random EUT
%     R_eut=.5; %(m)
%     n=5; %number of dipoles
%     
%     x=R_eut*rand(n,1)-R_eut/2;
%     y=R_eut*rand(n,1)-R_eut/2;
%     z=R_eut*rand(n,1)-R_eut/2;
%     tilt=pi*rand(n,1);
%     azimut=2*pi*rand(n,1);
%     
%     amplitude=rand(n,1);
%     phas=2*pi*rand(n,1);

    n=20;
    R_eut=.3; %m
    
    theta_eut=acos(2*rand(n,1)-1);
    phi_eut=2*pi*rand(n,1);
    x=R_eut*cos(phi_eut).*sin(theta_eut);
    y=R_eut*sin(phi_eut).*sin(theta_eut);
    z=R_eut*cos(theta_eut);
    tilt=pi*rand(n,1);
    azimut=2*pi*rand(n,1);
    
    amplitude=ones(n,1);%rand(n,1);
    phas=2*pi*rand(n,1);

    I = [x y z tilt azimut amplitude phas];
    
    
    %Free space caracterisation (perfect anechoic chamber)
    PHIac=2*pi*rand(M,1);
     THETAac=acos(2*rand(M,1)-1);%pi*rand(M,1);
    
    
    
    %for lll=1%:.2:6
    
    
    I1=I;
    
    
    for i=1:length(THETAac)
        
        X=R*cos(PHIac(i))*sin(THETAac(i));
        Y=R*sin(PHIac(i))*sin(THETAac(i));
        Z=R*cos(THETAac(i));
        DX = X-I1(:,1);
        DY = Y-I1(:,2);
        DZ = Z-I1(:,3);
        
        dist = sqrt(DX.^2+DY.^2+DZ.^2);
        
        dp=repmat(dist,1,length(f));
        fp=repmat(f,length(dist),1);
        phaseI=repmat(I1(:,7),1,length(f));
        phase=2*pi*dp.*fp/c+phaseI;
        ca    = cos(I1(:,4));
        sa    = sin(I1(:,4));
        cb    = cos(I1(:,5));
        sb    = sin(I1(:,5));
        
        clear fp
        distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
        disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
        distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;
        DXY=sqrt(DX.^2+DY.^2);
        clear DX DY DZ DXY
        distxy = sqrt(distx.^2+disty.^2);
        costheta = distz./dist;
        sintheta = distxy./dist;
        cosphi   = distx./distxy;
        sinphi   = disty./distxy;
        L =repmat(I1(:,6),1,length(f))*1./dp; %Amplitude & free space attenuation
        clear distx disty distz distxy
        clear dist dp
        %Projection in the usual rectangular coordinates
        Exx = sum(exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        Eyy = sum(exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        Ezz = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
        
        Ethac(i,o,:)= Exx*cos(THETAac(i))*cos(PHIac(i))+Eyy*cos(THETAac(i))*sin(PHIac(i))-Ezz*sin(THETAac(i));
        Ephac(i,o,:)= -Exx*sin(PHIac(i))+Eyy*cos(PHIac(i));
        Erac(i,o,:) = Exx*sin(THETAac(i))*cos(PHIac(i))+Eyy*sin(THETAac(i))*sin(PHIac(i))+Ezz*cos(THETAac(i));
        
        
        clear Exx Eyy Ezz
        
    end
    
    
    
    
    %%%AC with GP
    
    
    I2=I1;
    I2(:,3)=I1(:,3)+h+R_eut;
    Ip=I2;
    Ip(:,3)=-Ip(:,3);
    Ip(:,6)=Reflec*Ip(:,6);
    
    POS=[I2;Ip];

    
    for i=1:length(THETAac)
        
        X=R*cos(PHIac(i))*sin(THETAac(i));
        Y=R*sin(PHIac(i))*sin(THETAac(i));
        Z=R*cos(THETAac(i));
        %E-field calculation
        DX = X-POS(:,1);
        DY = Y-POS(:,2);
        DZ = Z-POS(:,3);
        
        dist = sqrt(DX.^2+DY.^2+DZ.^2);
        
        dp=repmat(dist,1,length(f));
        fp=repmat(f,length(dist),1);
        phaseI=repmat(POS(:,7),1,length(f));
        phase=2*pi*dp.*fp/c+phaseI;
        ca    = cos(POS(:,4));
        sa    = sin(POS(:,4));
        cb    = cos(POS(:,5));
        sb    = sin(POS(:,5));
        
        clear fp
        distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
        disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
        distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;
        DXY=sqrt(DX.^2+DY.^2);
        clear DX DY DZ DXY
        distxy = sqrt(distx.^2+disty.^2);
        costheta = distz./dist;
        sintheta = distxy./dist;
        cosphi   = distx./distxy;
        sinphi   = disty./distxy;
        L =repmat(POS(:,6),1,length(f))*1./dp; %Amplitude & free space attenuation
        clear distx disty distz distxy
        clear dist dp
        %Projection in the usual rectangular coordinates
        Exx = sum(exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        Eyy = sum(exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        Ezz = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
        
        Ethacgp(i,o,:)= Exx*cos(THETAac(i))*cos(PHIac(i))+Eyy*cos(THETAac(i))*sin(PHIac(i))-Ezz*sin(THETAac(i));
        Ephacgp(i,o,:)= -Exx*sin(PHIac(i))+Eyy*cos(PHIac(i));
        Eracgp(i,o,:) = Exx*sin(THETAac(i))*cos(PHIac(i))+Eyy*sin(THETAac(i))*sin(PHIac(i))+Ezz*cos(THETAac(i));
        
        clear Exx Eyy Ezz
        
    end
    
    
    
    %%%OATS
    
    Zr=Hmin+(Hmax-Hmin)*rand(M,1);
    %THETAoats=pi/2-atan(Zr./R);
    
    % [PHIac,THETAac]=meshgrid(phi,theta);
    
    for i=1:length(Zr)
        
        X=R*cos(PHIac(i));%*sin(THETAoats(i));
        Y=R*sin(PHIac(i));%*sin(THETAoats(i));
        Z=Zr(i);
        %E-field calculation
        DX = X-POS(:,1);
        DY = Y-POS(:,2);
        DZ = Z-POS(:,3);
        
        dist = sqrt(DX.^2+DY.^2+DZ.^2);
        
        dp=repmat(dist,1,length(f));
        fp=repmat(f,length(dist),1);
        phaseI=repmat(POS(:,7),1,length(f));
        phase=2*pi*dp.*fp/c+phaseI;
        ca    = cos(POS(:,4));
        sa    = sin(POS(:,4));
        cb    = cos(POS(:,5));
        sb    = sin(POS(:,5));
        
        clear fp
        distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
        disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
        distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;
        DXY=sqrt(DX.^2+DY.^2);
        clear DX DY DZ DXY
        distxy = sqrt(distx.^2+disty.^2);
        costheta = distz./dist;
        sintheta = distxy./dist;
        cosphi   = distx./distxy;
        sinphi   = disty./distxy;
        L =repmat(POS(:,6),1,length(f))*1./dp; %Amplitude & free space attenuation
        clear distx disty distz distxy
        clear dist dp
        %Projection in the usual rectangular coordinates
        Exx = sum(exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        Eyy = sum(exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
        %Ezz = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
       
        Ezoats(i,o,:) = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
        Ephoats(i,o,:) = -sin(PHIac(i)).*Exx+cos(PHIac(i)).*Eyy;
        Eroats(i,o,:) = cos(PHIac(i)).*Exx+sin(PHIac(i)).*Eyy;
        
        clear Exx Eyy Ezz
        
    end
    
    
    if mod(o,10)==0
        disp(o)
        toc
        
    end
    
    
end

Max_ac=max(squeeze(max([abs(Ethac);abs(Ephac)])));
Max_acgp=max(squeeze(max([abs(Ethacgp);abs(Ephacgp)])));
Max_oats=max(squeeze(max([abs(Ezoats);abs(Ephoats)])));

max_ac=mean(squeeze(max([abs(Ethac);abs(Ephac)])));
max_acgp=mean(squeeze(max([abs(Ethacgp);abs(Ephacgp)])));
max_oats=mean(squeeze(max([abs(Ezoats);abs(Ephoats)])));

mean_ac=mean(squeeze(mean([abs(Ethac);abs(Ephac)])));
mean_acgp=mean(squeeze(mean([abs(Ethacgp);abs(Ephacgp)])));
mean_oats=mean(squeeze(mean([abs(Ezoats);abs(Ephoats)])));



figure(1)
plot(f/1e6,mean_acgp./mean_ac,f/1e6,max_acgp./max_ac,f/1e6,mean_oats./mean_ac,f/1e6,max_oats./max_ac)%,f/1e6,Max_acgp./Max_ac,f/1e6,Max_oats./Max_ac)
legend('AC_{GP}m/ACm','AC_{GP}max/ACmax','OATSm/ACm','OATSmax/ACmax')
grid on