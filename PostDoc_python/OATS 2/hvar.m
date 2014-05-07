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
f=10e6:10e6:1e9;
lambda=c./f;

%EUT
%positions of the dipoles and amplitude and phase

%random EUT
% L=.5; %(m)
% n=5; %number of dipoles
% 
% x=L*rand(n,1)-L/2;
% y=L*rand(n,1)-L/2;
% z=L*rand(n,1)+h;
% tilt=pi*rand(n,1);
% azimut=2*pi*rand(n,1);
% 
% amplitude=rand(n,1);
% phas=2*pi*rand(n,1);

% x=[0;0]%[-.8;-.6;-.4;-.2;.2;.4;.6;.8]
% y=[0;0]%[-.8;-.6;-.4;-.2;.2;.4;.6;.8]
% z=[.3;-.3]+h;%[h;h;h;h;h;h;h;h]
% tilt=[0;0]%[0;0;0;0;0;0;0;0]
% azimut=[0;0]%[0;0;0;0;0;0;0;0]
% amplitude=[1;1];%[1;1;1;1;1;1;1;1]
% phas=[0;0]%[0;0;0;0;0;0;0;0]

% I = [x y z tilt azimut amplitude phas];

load antrand5_50cm.mat


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
        DX = X-I(:,1);
        DY = Y-I(:,2);
        DZ = Z-I(:,3); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dist = sqrt(DX.^2+DY.^2+DZ.^2);
        
        dp=repmat(dist,1,length(f));
        fp=repmat(f,length(dist),1);
        phaseI=repmat(I(:,7),1,length(f));
        phase=2*pi*dp.*fp/c+phaseI;
        ca    = cos(I(:,4));
        sa    = sin(I(:,4));
        cb    = cos(I(:,5));
        sb    = sin(I(:,5));
        
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
        L =repmat(I(:,6),1,length(f))*1./dp; %Amplitude & free space attenuation
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
    Dac(i) = 4*pi / omega;
end

u=0;
for h=.1:.1:2
    load antrand5_50cm.mat
u=u+1;
H(u)=h
%%%OATS
I(:,3)=I(:,3)+h;
Ip = I;
Ip(:,3)=-Ip(:,3);
Ip(:,6)=Reflec*Ip(:,6);

POS=[I;Ip];

tic
%Measurement positions
%phi=(1:5:360)*pi/180;
Zr=Hmin:.05:Hmax;

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
    Doats(u,i,:) = 4*pi / omega;
end
end


save('result_n5_50cm-H0-2m.mat','H','f','Dac','Doats')

for i=1:20
    rD(i,:)=Doats(i,:)./Dac;
end
for i=1:20
plot(f,rD(i,:))
ylim([0 70])
getframe;

end
% figure(2)
% for i=1:length(f)
%     A=Ethac(:,:,i);
%     hist(abs(A(:)),50)
%     title(num2str(f(i)/1e6))
%     getframe;
% end
% 
% 
% 
% figure(3)
% maxca=max(max(max(max([abs(Ethac);abs(Ephac);abs(Erac)]))))
% cmin=0;
% for i=1:length(f)
%     subplot(1,3,1)
%     imagesc(theta,phi,abs(Ethac(:,:,i)'))
%     title(['E_\theta [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     caxis([0 maxca])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     xlabel('\theta [rad]')
%     ylabel('\phi [rad]')
%     
%     subplot(1,3,2)
%     imagesc(theta,phi,abs(Ephac(:,:,i)'))
%     title(['E_\phi [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     
%     caxis([0 maxca])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     xlabel('\theta [rad]')
%     ylabel('\phi [rad]')
%     
%     subplot(1,3,3)
%     imagesc(theta,phi,abs(Erac(:,:,i)'))
%     title(['E_r [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     caxis([0 maxca])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     xlabel('\theta [rad]')
%     ylabel('\phi [rad]')
%     getframe;
% end
% 
% figure(4)
% maxoats=max(max(max(max([abs(Ez);abs(Ephi);abs(Er)]))))
% cmin=0;
% for i=1:length(f)
%     subplot(1,3,1)
%     imagesc(Zr,phi,abs(Ez(:,:,i)'))
%     title(['E_z [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     caxis([0 maxoats])
%     xlabel('z [m]')
%     ylabel('\phi [rad]')
%     
%     subplot(1,3,2)
%     imagesc(Zr,phi,abs(Ephi(:,:,i)'))
%     title(['E_\phi [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     caxis([0 maxoats])
%     xlabel('z [m]')
%     ylabel('\phi [rad]')
%     
%     subplot(1,3,3)
%     imagesc(Zr,phi,abs(Er(:,:,i)'))
%     title(['E_r [V/m], f=',num2str(f(i)/1e6),' MHz'])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     caxis([0 maxoats])
%     xlabel('z [m]')
%     ylabel('\phi [rad]')
%     
%     getframe;
% end
% 
% figure(5)
% for i=1:length(f)
%     
%     imagesc(phi,Zr,abs(Ez(:,:,i).^2+Ephi(:,:,i).^2))
%     title(['P [mW], f=',num2str(f(i)/1e6),' MHz'])
%     colorbar('location','southoutside')
%     axis equal
%     axis tight
%     xlabel('\phi [rad]')
%     ylabel('z [m]')
%     getframe;
% end
% 
% 



