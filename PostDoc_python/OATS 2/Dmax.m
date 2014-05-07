clear java
clear all
c=2.998e8;

%open area test setup
Hmax = 4;%maximum heigth of the antenna (m)
Hmin = 1;%minimum heigth of the antenna (m)
R = 10;  %radius (m)
h = 1;   %heigth of the EUT (m)
Reflec = 1; %Reflection coefficient
%signal
f=1e6:1e6:3e6%:100e6:1e9;
lambda=c./f;


%M=50;
MC=30;
%EUT
%positions of the dipoles and amplitude and phase



%Free space caracterisation (perfect anechoic chamber)
phi=(1:1:360)*pi/180;
theta=(1:1:180)*pi/180;%acos(-1:2/180:1);

[PHIac,THETAac]=meshgrid(phi,theta);
tic
n=1;
%for u=1:length(RR)
%EUT
%positions of the dipoles and amplitude and phase
tic
for o=1:MC
    
    R_eut=0.3;%RR(u); %m
    
    theta_eut=acos(2*rand(n,1)-1);%pi*rand(n,1);
    phi_eut=2*pi*rand(n,1);
    x=R_eut*cos(phi_eut).*sin(theta_eut);
    y=R_eut*sin(phi_eut).*sin(theta_eut);
    z=R_eut*cos(theta_eut);
    tilt=pi*rand(n,1);
    azimut=2*pi*rand(n,1);
    ld=.01;%dipole length
    
    amplitude=ones(n,1);%rand(n,1);
    phas=2*pi*rand(n,1);
    
    I = [x y z tilt azimut amplitude phas];
    
    for i=1:length(theta)
        for j=1:1:length(phi)
            X=R*cos(PHIac(i,j))*sin(THETAac(i,j));
            Y=R*sin(PHIac(i,j))*sin(THETAac(i,j));
            Z=R*cos(THETAac(i,j));
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
            L =repmat(I(:,6),1,length(f))*1./dp;%*ld*c./fp/2*377; %Amplitude & free space attenuation
            clear distx disty distz distxy
            clear dist dp
            
            %Projection in the usual rectangular coordinates
            Exx = (exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
            Eyy = (exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
            Ezz = (exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
            
            Ethac(i,j,:)= Exx*cos(THETAac(i))*cos(PHIac(i))+Eyy*cos(THETAac(i))*sin(PHIac(i))-Ezz*sin(THETAac(i));
            Ephac(i,j,:)= -Exx*sin(PHIac(i))+Eyy*cos(PHIac(i));
            Erac(i,j,:) = Exx*sin(THETAac(i))*cos(PHIac(i))+Eyy*sin(THETAac(i))*sin(PHIac(i))+Ezz*cos(THETAac(i));
        end
        
        
        
        
        
    end
    nt = length(THETAac(1,:));
    np = length(PHIac(:,1));
    
    dtheta = pi/nt;
    dphi = (2*pi)/np;
    
    for i=1:length(f)
        Fa2=abs(Ephac(:,:,i)).^2+abs(Ethac(:,:,i)).^2;
        Fa=Fa2/max(max(Fa2));
        omega = sum(sum(Fa .* sin(THETAac)*dtheta*dphi));
        Dac(o,i) = 4*pi / omega;
    end
    clear Exx Eyy Ezz Ephac Ethac Erac
    
    %if mod(o,10)==0
    disp(o)
    toc
    
    %end
    
end



ka=2*pi./lambda*R_eut;

k=1:.1:6;
Dmaxth=1/2*(0.577+log((4*k).^2+8*k)+1./(8*k.^2+16*k));
kk=5:.1:10;
Dmaxth2=0.982+log(kk);


plot(ka,mean(Dac),k,Dmaxth)



