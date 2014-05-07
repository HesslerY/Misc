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
f=10e6:20e6:3e9;
lambda=c./f;

%EUT
%positions of the dipoles and amplitude and phase

%random EUT
D=.5; %(m)
n=5; %number of dipoles

x=D*rand(n,1)-D/2;
y=D*rand(n,1)-D/2;
z=D*rand(n,1)-D/2;
tilt=pi*rand(n,1);
azimut=2*pi*rand(n,1);

amplitude=rand(n,1);
phas=2*pi*rand(n,1);

% x=[0;0]%[-.8;-.6;-.4;-.2;.2;.4;.6;.8]
% y=[0;0]%[-.8;-.6;-.4;-.2;.2;.4;.6;.8]
% z=[.3;-.3]+h;%[h;h;h;h;h;h;h;h]
% tilt=[0;0]%[0;0;0;0;0;0;0;0]
% azimut=[0;0]%[0;0;0;0;0;0;0;0]
% amplitude=[1;1];%[1;1;1;1;1;1;1;1]
% phas=[0;0]%[0;0;0;0;0;0;0;0]

I = [x y z tilt azimut amplitude phas];
k=0
for lll=1:.1:6
    k=k+1;
    disp(k/length(1:.05:6));
    DD(k)=D*lll;
    I1=I;
    I1(:,1)=I1(:,1)*lll;
    I1(:,2)=I1(:,2)*lll;
    I1(:,3)=I1(:,3)*lll;
    
    
    %load antrand5_50cm.mat
    
    
    %Free space caracterisation (perfect anechoic chamber)
    phi=(1:2:360)*pi/180;
    theta=(1:2:180)*pi/180;
    
    [PHIac,THETAac]=meshgrid(phi,theta);
    tic
    for i=1:length(theta)
        for j=1:length(phi)
            X=R*cos(PHIac(i,j))*sin(THETAac(i,j));
            Y=R*sin(PHIac(i,j))*sin(THETAac(i,j));
            Z=R*cos(THETAac(i,j));
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
            
            Ethac(i,j,:)= Exx*cos(THETAac(i,j))*cos(PHIac(i,j))+Eyy*cos(THETAac(i,j))*sin(PHIac(i,j))-Ezz*sin(THETAac(i,j));
            Ephac(i,j,:)= -Exx*sin(PHIac(i,j))+Eyy*cos(PHIac(i,j));
            Erac(i,j,:) = Exx*sin(THETAac(i,j))*cos(PHIac(i,j))+Eyy*sin(THETAac(i,j))*sin(PHIac(i,j))+Ezz*cos(THETAac(i,j));
            
            
            clear Exx Eyy Ezz
            
        end
    end
    toc
    
    nt = length(THETAac(1,:));
    np = length(PHIac(:,1));
    
    dtheta = pi/nt;
    dphi = (2*pi)/np;
    
    for i=1:length(f)
        Fa2=abs(Ephac(:,:,i)).^2+abs(Ethac(:,:,i)).^2;
        Fa=Fa2/max(max(Fa2));
        omega = sum(sum((Fa) .* sin(THETAac)*dtheta*dphi));
        Dac(i,k,:) = 4*pi / omega;
    end
    
    %%%OATS
    %%%OATS
    %%%OATS
    %%%OATS
    
    
    I2=I1;
    I2(:,3)=I1(:,3)+h+D/2;
    Ip=I2;
    Ip(:,3)=-Ip(:,3);
    Ip(:,6)=Reflec*Ip(:,6);
    
    POS=[I2;Ip];
    
    tic
    %Measurement positions
    %phi=(1:5:360)*pi/180;
    Zr=Hmin:.1:Hmax;
    
    [PHI,Z]=meshgrid(phi,Zr);
    tic
    for i=1:length(Zr)
        for j=1:length(phi)
            X=R*cos(PHI(i,j));
            Y=R*sin(PHI(i,j));
            %E-field calculation
            DX = X-POS(:,1);
            DY = Y-POS(:,2);
            DZ = Z(i,j)-POS(:,3);
            
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
            
            Ez(i,j,:) = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
            Ephi(i,j,:) = -sin(PHI(i,j)).*Exx+cos(PHI(i,j)).*Eyy;
            Er(i,j,:) = cos(PHI(i,j)).*Exx+sin(PHI(i,j)).*Eyy;
            Ex(i,j,:) = Exx;
            Ey(i,j,:) = Eyy;
            
            clear Exx Eyy
            
        end
    end
    toc
    
    
    nt = length(Z(1,:));
    np = length(PHI(:,1));
    
    dtheta = pi/nt;
    dphi = (2*pi)/np;
    
    for i=1:length(f)
        Fa2=abs(Ephi(:,:,i)).^2+abs(Er(:,:,i)).^2+abs(Ez(:,:,i)).^2;
        Fa=Fa2/max(max(Fa2));
        omega = sum(sum((Fa) .* sin(Z)*dtheta*dphi));
        Doats(i,k,:) = 4*pi / omega;
    end
    
end


plot(f,Dac)%,f,Doats)

for i=1:51
    plot(f,Doats(:,i))
    ylim([0 25])
    getframe
    pause(.2)
end

