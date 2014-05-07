clear java
clear all
c=2.998e8;

M=50;
MC=100;
%open area test setup
Hmax = 4;%maximum heigth of the antenna (m)
Hmin = 1;%minimum heigth of the antenna (m)
R = 10;  %radius (m)
h = 1;   %heigth of the EUT (m)
Reflec = 1; %Reflection coefficient
%signal
f=1.e8:1e8:10e9;
lambda=c./f;

% Ephac=zeros(M,MC),length(f));
% Ephacgp=zeros(M,MC),length(f));
% Ephoats=zeros(M,MC),length(f));
%
% Erac=zeros(M,MC,length(f));
% Eracgp=zeros(M,MC,length(f));
% Eroats=zeros(M,MC,length(f));
% Ethac=zeros(M,MC,length(f));
% Ethacgp=zeros(M,MC,length(f));
% Ezoats=zeros(M,MC,length(f));
%RR=[.01;0.05;0.1;0.3;0.5;1];
n=5;
%for u=1:length(RR)
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
        
        %n=20;
        R_eut=0.3%RR(u); %m
        
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
            
            Ethac(u,i,o,:)= Exx*cos(THETAac(i))*cos(PHIac(i))+Eyy*cos(THETAac(i))*sin(PHIac(i))-Ezz*sin(THETAac(i));
            Ephac(u,i,o,:)= -Exx*sin(PHIac(i))+Eyy*cos(PHIac(i));
            Erac(u,i,o,:) = Exx*sin(THETAac(i))*cos(PHIac(i))+Eyy*sin(THETAac(i))*sin(PHIac(i))+Ezz*cos(THETAac(i));
            
            
            clear Exx Eyy Ezz
            
        end
        
        
        
            if mod(o,10)==0
                disp(o)
                toc
        
            end
        
        
%    % end
%     disp(o)
%     toc
end


% for i=1:length(RR)
%     for j=1:MC
%         mean_ac(i,j)=mean(mean(squeeze([abs(Ethac(i,:,j));abs(Ephac(i,:,j))])));
%         
%     end
% end
% 
% m975=prctile(mean_ac,97.5,2);
% m25=prctile(mean_ac,2.5,2);
% 
% 
% plot(RR,mean(mean_ac'),RR,m975,RR,m25)
% legend('$<E_R(N)>$ (6)','$<E_R(N)>$ MC simulations','$95~\%$ CI of $<E_R(N)>$')
% xlabel('$R$')
% ylabel('[V/m]')
% 
% plot(1:20,mean_th',1:20,mean_ac)
% %plot(1:20,mean_th'/max(mean_th),1:20,mean_ac/max(mean_ac))
% 
% 
% %plot(1:20,1./(mean_ac./mean_th))

X=0:.001:4
th=cdf('rayleigh',X,1)

%
%%%STAT
figure(2)
hold on
plot(X,th,'--r','LineWidth',2)
for i=1
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    
end
xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=1$ cm, 50 simulations','Interpreter','Latex','Fontsize',12)

h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)
grid on
ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)
grid on



figure(3)
hold on
plot(X,th,'--r','LineWidth',2)
for i=2
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    
    
end
xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=5$ cm,, 50 simulations','Interpreter','Latex','Fontsize',12)
grid on
h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)
grid on
ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)
grid on


figure(4)
hold on
plot(X,th,'--r','LineWidth',2)
for i=3
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    
    
end
xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=10$ cm, 50 simulations','Interpreter','Latex','Fontsize',12)
grid on
h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)
grid on
ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)
grid on


figure(5)
hold on
plot(X,th,'--r','LineWidth',2)
for i=4
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    
    
end
xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=30$ cm, 50 simulations','Interpreter','Latex','Fontsize',12)
grid on
h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)
grid on
ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)
grid on



figure(6)
hold on
plot(X,th,'--r','LineWidth',2)
for i=5
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    
    
end

xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=50$ cm, 50 simulations','Interpreter','Latex','Fontsize',12)
grid on
h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)
grid on
ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)
grid on


figure(7)
hold on
plot(X,th,'--r','LineWidth',2)
for i=6
    for j=1:1:MC
        A=[abs(Ethac(i,:,j));abs(Ephac(i,:,j))];
        B=A(:)/mean(A(:));
        [f,x] = ecdf(B/mean(B)*sqrt(pi/2));
        stairs(x,f)
        
    end
    plot(X,th,'--r','LineWidth',2)
    title('$N=20$, 50 simulations')
    
end
xlim([0 4])
plot(X,th,'--r','LineWidth',2)
title('$R=1$ m, 50 simulations','Interpreter','Latex','Fontsize',12)
grid on
h=legend('Theoretical Rayleigh CDF','MC experiments');
set(h,'Interpreter','Latex')
xlabel('$\frac{E_R}{<E_R>}\sqrt{\pi/2}$','Interpreter','Latex','Fontsize',12)

ylabel('Cumulative probability','Interpreter','Latex','Fontsize',12)

