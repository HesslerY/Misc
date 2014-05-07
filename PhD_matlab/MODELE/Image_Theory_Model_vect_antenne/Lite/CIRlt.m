%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                               %%%%%%
%%%%%            CHANNEL IMPULSE RESPONSE FUNCTION Lite             %%%%%%
%%%%%                                                               %%%%%%
%%%%%           by E. Amador (emmanuel.amador@gmail.com)            %%%%%%
%%%%%                         IETR/DGA                              %%%%%%
%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Sx,Sy,Sz]=CIRlt(X_1,Y_1,Z_1)
d=1; %value of the elementary pulse

global Lt c Rx Ry Rz N POS
%Channel Impulse Response
Sx=zeros(1,N);
Sy=zeros(1,N);
Sz=zeros(1,N);
for j=1:1:length(POS) %loop over the image-currents... to be vctorized... one day may be...
    
    DX=X_1-POS(j,1);
    DY=Y_1-POS(j,2);
    DZ=Z_1-POS(j,3);
    
    dist=sqrt(DX^2+DY^2+DZ^2);
    zl=round((N-1)*dist/c/Lt);
    
    if zl<N
        alpha=POS(j,5);
        beta=POS(j,6);
        ca=cos(alpha);
        sa=sin(alpha);
        cb=cos(beta);
        sb=sin(beta);
        
        %Ralpha/beta, coordinates transformation matrix:
        %                           Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
        %                                     -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
        %                                     -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];
        
        %rectangular coordinates calculation in the local system attached to the considered current (developped expressions, a matrix product is slower)
        distx=((-sb)^2+(1-(-sb)^2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ;
        disty=(-sb*cb*(1-ca))*DX+((cb)^2+(1-cb^2)*ca)*DY+(sb*sa)*DZ;
        distz=(-cb*sa)*DX+(-sb*sa)*DY+(ca)*DZ;
        
        distxy=sqrt(distx^2+disty^2);
        
        R=Rx^POS(j,4)%Attenuation
        E=d/dist; % Free-space attenuation
        costheta=distz/dist;
        sintheta=distxy/dist;
        cosphi=distx/distxy;
        sinphi=disty/distxy;
        Antth=-sintheta; %dipole radiation pattern
        
        %Reverse transformation matrix
        %                         Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
        %                                     (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
        %                                     (cb*sa) (sb*sa) ca];
        
      
        %Projection in the usual rectangular coordinates
        
        Vx=(((-sb)^2+(1-(-sb)^2)*ca)*(Antth*costheta*cosphi)+(-sb*cb*(1-ca))*(Antth*costheta*sinphi)+(-cb*sa)*(-sintheta*Antth));
        Vy=((-sb*cb*(1-ca))*(Antth*costheta*cosphi)+((cb)^2+(1-(cb)^2)*ca)*(Antth*costheta*sinphi)+(-sb*sa)*(-sintheta*Antth));
        Vz=((cb*sa)*(Antth*costheta*cosphi)+(sb*sa)*(Antth*costheta*sinphi)+ca*(-sintheta*Antth));
        
        
        %Three axis channel impulse response construction
        Sx(zl+1)=Sx(zl+1)+R*E*Vx;
        Sy(zl+1)=Sy(zl+1)+R*E*Vy;
        Sz(zl+1)=Sz(zl+1)+R*E*Vz;
        
    end
end
