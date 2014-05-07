% clear java
% clear all
% c=2.998e8;
% 
% M=50;
% MC=500;
% %open area test setup
% Hmax = 4;%maximum heigth of the antenna (m)
% Hmin = 1;%minimum heigth of the antenna (m)
% R = 10;  %radius (m)
% h = 1;   %heigth of the EUT (m)
% Reflec = 1; %Reflection coefficient
% %signal
% f=1e9;
% lambda=c./f;
% 
% % Ephac=zeros(M,MC),length(f));
% % Ephacgp=zeros(M,MC),length(f));
% % Ephoats=zeros(M,MC),length(f));
% %  
% % Erac=zeros(M,MC,length(f));
% % Eracgp=zeros(M,MC,length(f));
% % Eroats=zeros(M,MC,length(f));
% % Ethac=zeros(M,MC,length(f));
% % Ethacgp=zeros(M,MC,length(f));
% % Ezoats=zeros(M,MC,length(f));
% 
% for n=1:20;
% %EUT
% %positions of the dipoles and amplitude and phase
% tic
% for o=1:MC
%     %random EUT
% %     R_eut=.5; %(m)
% %     n=5; %number of dipoles
% %     
% %     x=R_eut*rand(n,1)-R_eut/2;
% %     y=R_eut*rand(n,1)-R_eut/2;
% %     z=R_eut*rand(n,1)-R_eut/2;
% %     tilt=pi*rand(n,1);
% %     azimut=2*pi*rand(n,1);
% %     
% %     amplitude=rand(n,1);
% %     phas=2*pi*rand(n,1);
% 
%     %n=20;
%     R_eut=.3; %m
%     
%     theta_eut=pi*rand(n,1);
%     phi_eut=2*pi*rand(n,1);
%     x=R_eut*cos(phi_eut).*sin(theta_eut);
%     y=R_eut*sin(phi_eut).*sin(theta_eut);
%     z=R_eut*cos(theta_eut);
%     tilt=pi*rand(n,1);
%     azimut=2*pi*rand(n,1);
%     
%     amplitude=ones(n,1);%rand(n,1);
%     phas=2*pi*rand(n,1);
% 
%     I = [x y z tilt azimut amplitude phas];
%     
%     
%     %Free space caracterisation (perfect anechoic chamber)
%     PHIac=2*pi*rand(M,1);
%     THETAac=pi*rand(M,1);
%     
%     
%     
%     %for lll=1%:.2:6
%     
%     
%     I1=I;
%     
%     
%     for i=1:length(THETAac)
%         
%         X=R*cos(PHIac(i))*sin(THETAac(i));
%         Y=R*sin(PHIac(i))*sin(THETAac(i));
%         Z=R*cos(THETAac(i));
%         DX = X-I1(:,1);
%         DY = Y-I1(:,2);
%         DZ = Z-I1(:,3);
%         
%         dist = sqrt(DX.^2+DY.^2+DZ.^2);
%         
%         dp=repmat(dist,1,length(f));
%         fp=repmat(f,length(dist),1);
%         phaseI=repmat(I1(:,7),1,length(f));
%         phase=2*pi*dp.*fp/c+phaseI;
%         ca    = cos(I1(:,4));
%         sa    = sin(I1(:,4));
%         cb    = cos(I1(:,5));
%         sb    = sin(I1(:,5));
%         
%         clear fp
%         distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
%         disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
%         distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;
%         DXY=sqrt(DX.^2+DY.^2);
%         clear DX DY DZ DXY
%         distxy = sqrt(distx.^2+disty.^2);
%         costheta = distz./dist;
%         sintheta = distxy./dist;
%         cosphi   = distx./distxy;
%         sinphi   = disty./distxy;
%         L =repmat(I1(:,6),1,length(f))*1./dp; %Amplitude & free space attenuation
%         clear distx disty distz distxy
%         clear dist dp
%         %Projection in the usual rectangular coordinates
%         Exx = sum(exp(1i*phase).*L.*repmat(((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
%         Eyy = sum(exp(1i*phase).*L.*repmat((((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta)))),1,length(f)));
%         Ezz = sum(exp(1i*phase).*L.*repmat((((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta)))),1,length(f)));
%         
%         Ethac(n,i,o,:)= Exx*cos(THETAac(i))*cos(PHIac(i))+Eyy*cos(THETAac(i))*sin(PHIac(i))-Ezz*sin(THETAac(i));
%         Ephac(n,i,o,:)= -Exx*sin(PHIac(i))+Eyy*cos(PHIac(i));
%         Erac(n,i,o,:) = Exx*sin(THETAac(i))*cos(PHIac(i))+Eyy*sin(THETAac(i))*sin(PHIac(i))+Ezz*cos(THETAac(i));
%         
%         
%         clear Exx Eyy Ezz
%         
%     end
%     
% 
%     
% %     if mod(o,10)==0
% %         disp(o)
% %         toc
% %         
% %     end
%     
%     
% end
% disp(n)
% toc
% end


for i=1:20
mean_th(i)=1/R*sqrt(i*pi/2)
%Max_ac(i)=max(squeeze(max([abs(Ethac);abs(Ephac)])));
%max_ac=mean(squeeze(max([abs(Ethac);abs(Ephac)])));
mean_ac(i)=mean(mean(mean(squeeze([abs(Ethac(i,:,:));abs(Ephac(i,:,:))]))));
end
%plot(1:20,mean_th',1:20,mean_ac)
plot(1:20,mean_th'/max(mean_th),1:20,mean_ac/max(mean_ac))


plot(1:20,1./(mean_ac./mean_th))
