clear java
clear all
c=299792458;

R = 100;  %radius (m)
f=1e9;%1e6:5e6:3e9;
lambda=c./f;
h=1;
Reflec=1;
M=30;
MC=500;
%EUT


%[PHIac,THETAac]=meshgrid(phi,theta);
R_eut=.3; %m
tic
%random EUT
for n=1:20; %number of dipoles
    tic
    for o=1:MC
        phi=2*pi*rand(M,1);
        theta=acos(2*rand(M,1)-1);
        
        
        theta_eut=acos(2*rand(n,1)-1);
        phi_eut=2*pi*rand(n,1);
        x=R_eut*cos(phi_eut).*sin(theta_eut);
        y=R_eut*sin(phi_eut).*sin(theta_eut);
        z=R_eut*cos(theta_eut);
        tilt=acos(2*rand(n,1)-1);%pi/2*ones(n,1);
        azimut=2*pi*rand(n,1);
        
        ld=.1;
        
        amplitude=ones(n,1);
        phas=2*pi*rand(n,1);
        
        I = [x y z tilt azimut amplitude phas];
        
        
        %Free space caracterisation (perfect anechoic chamber)
        
        if n==1
            for i=1:length(theta)
                
                X=R*cos(phi(i))*sin(theta(i));
                Y=R*sin(phi(i))*sin(theta(i));
                Z=R*cos(theta(i));
                DX = X-I(:,1);
                DY = Y-I(:,2);
                DZ = Z-I(:,3);
                
                dist = sqrt(DX.^2+DY.^2+DZ.^2);
                
                dp=repmat(dist,1,length(f));
                fp=repmat(f,length(dist),1);
                phaseI=repmat(I(:,7),1,length(f));
                phase=2*pi*dp.*fp/c+phaseI;
                ca    = cos(I(:,4));
                sa    = sin(I(:,4));
                cb    = cos(I(:,5));
                sb    = sin(I(:,5));
                
                
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
                L =repmat(I(:,6),1,length(f))*1./dp*ld.*fp/c/2*377; %Amplitude & free space attenuation
                clear distx disty distz distxy
                clear dist dp
                %Projection in the usual rectangular coordinates
                Exx = (exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
                Eyy = (exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
                Ezz = (exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
                
                Ethac(n,o,i,:)= Exx*cos(theta(i))*cos(phi(i))+Eyy*cos(theta(i))*sin(phi(i))-Ezz*sin(theta(i));
                Ephac(n,o,i,:)= -Exx*sin(phi(i))+Eyy*cos(phi(i));
                %Erac(n,o,i,:) = Exx*sin(theta(i))*cos(phi(i))+Eyy*sin(theta(i))*sin(phi(i))+Ezz*cos(theta(i));
                
            end
            clear Exx Eyy Ezz
        else
            for i=1:length(theta)
                
                X=R*cos(phi(i))*sin(theta(i));
                Y=R*sin(phi(i))*sin(theta(i));
                Z=R*cos(theta(i));
                DX = X-I(:,1);
                DY = Y-I(:,2);
                DZ = Z-I(:,3);
                
                dist = sqrt(DX.^2+DY.^2+DZ.^2);
                
                dp=repmat(dist,1,length(f));
                fp=repmat(f,length(dist),1);
                phaseI=repmat(I(:,7),1,length(f));
                phase=2*pi*dp.*fp/c+phaseI;
                ca    = cos(I(:,4));
                sa    = sin(I(:,4));
                cb    = cos(I(:,5));
                sb    = sin(I(:,5));
                
                
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
                L =repmat(I(:,6),1,length(f))*1./dp*ld.*fp/c/2*377;%377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld./dp; %Amplitude & free space attenuation
                clear distx disty distz distxy
                clear dist dp
                %Projection in the usual rectangular coordinates
                Exx = sum(exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
                Eyy = sum(exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
                Ezz = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
                
                Ethac(n,o,i,:)= Exx*cos(theta(i))*cos(phi(i))+Eyy*cos(theta(i))*sin(phi(i))-Ezz*sin(theta(i));
                Ephac(n,o,i,:)= -Exx*sin(phi(i))+Eyy*cos(phi(i));
                %Erac(n,o,i,:) = Exx*sin(theta(i))*cos(phi(i))+Eyy*sin(theta(i))*sin(phi(i))+Ezz*cos(theta(i));
                
                
                clear Exx Eyy Ezz
                
            end
        end
    end
    toc
end
toc

for i=1:20
    A=[Ethac(i,:,:);Ethac(i,:,:)];
    E(i,:)=mean(mean(abs(A),3),1);
end


n=1:20;

Erth=377*2*pi*f/c*ld/16/R*sqrt(n*pi/2)/2;

plot(n,mean(E'),n,Erth,n,prctile(E',2.5),n,prctile(E',97.5))
