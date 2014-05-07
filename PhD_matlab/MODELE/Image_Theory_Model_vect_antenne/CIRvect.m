%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                               %%%%%%
%%%%%              CHANNEL IMPULSE RESPONSE FUNCTION                %%%%%%
%%%%%                                                               %%%%%%
%%%%%           by E. Amador (emmanuel.amador@gmail.com)            %%%%%%
%%%%%                         IETR/DGA                              %%%%%%
%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Sx,Sy,Sz]=CIRvect(X_1,Y_1,Z_1)
d=1; %value of the elementary pulse
tic
global Lt c Rx Ry Rz N POS
%Channel Impulse Response
Sx=zeros(1,N);
Sy=zeros(1,N);
Sz=zeros(1,N);
%for j=1:1:length(POS) %loop over the image-currents... to be vctorized... one day may be...
    
    DX=X_1-POS(:,1);
    DY=Y_1-POS(:,2);
    DZ=Z_1-POS(:,3);
    
    dist=sqrt(DX.^2+DY.^2+DZ.^2);
    zl=dist./c;
    zn=round((N-1)*dist./c/Lt);
    
    
        alpha=POS(:,9);
        beta=POS(:,10);
        ca=cos(alpha);
        sa=sin(alpha);
        cb=cos(beta);
        sb=sin(beta);
        
        %Ralpha/beta, coordinates transformation matrix:
        %                           Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
        %                                     -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
        %                                     -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];
        
        %rectangular coordinates calculation in the local system attached to the considered current (developped expressions, a matrix product is slower)
        distx=((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
        disty=(-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
        distz=(-cb.*sa).*DX+(-sb.*sa).*DY+(ca).*DZ;
        
        distxy=sqrt(distx.^2+disty.^2);
        
       
        E=Rx.^POS(:,4).*Ry.^POS(:,5).*Rz.^POS(:,6).*POS(:,8).*d./dist; % Free-space attenuation
        costheta=distz./dist;
        sintheta=distxy./dist;
        cosphi=distx./distxy;
        sinphi=disty./distxy;
        Antth=-sintheta; %dipole radiation pattern
        
        %Reverse transformation matrix
        %                         Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
        %                                     (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
        %                                     (cb*sa) (sb*sa) ca];
        
      
    
        
        %Three axis channel impulse response construction
        Sxb=E.*(((-sb).^2+(1-(-sb).^2).*ca).*(Antth.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(Antth.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*Antth));
        Syb=E.*((-sb.*cb.*(1-ca)).*(Antth.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(Antth.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*Antth));
        Szb=E.*((cb.*sa).*(Antth.*costheta.*cosphi)+(sb.*sa).*(Antth.*costheta.*sinphi)+ca.*(-sintheta.*Antth));

toc
Sx=zeros(1,N);
Sy=zeros(1,N);
Sz=zeros(1,N);
for j=1:1:length(Sxb)
    if zn(j)<N+1
Sx(zn(j))=Sx(zn(j))+Sxb(j);
Sy(zn(j))=Sy(zn(j))+Syb(j);
Sz(zn(j))=Sz(zn(j))+Szb(j);
    end
end
toc
