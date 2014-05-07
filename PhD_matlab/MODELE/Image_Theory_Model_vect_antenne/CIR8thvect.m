
%CIR8th.m computes the CIR of a given 8th of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%        CHANNEL IMPULSE RESPONSE Eighth FUNCTION V 2.0         %
%                                                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sx8th,Sy8th,Sz8th,dist,azim,elev]=CIR8thvect(X_1,Y_1,Z_1)

global Lt c R N POS


DX = X_1-POS(:,1);
DY = Y_1-POS(:,2);
DZ = Z_1-POS(:,3);

dist = sqrt(DX.^2+DY.^2+DZ.^2);

ca    = cos(POS(:,5));
sa    = sin(POS(:,5));
cb    = cos(POS(:,6));
sb    = sin(POS(:,6));

distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;

distxy = sqrt(distx.^2+disty.^2);

L = R.^POS(:,4); %Attenuation

costheta = distz./dist;
sintheta = distxy./dist;
cosphi   = distx./distxy;
sinphi   = disty./distxy;
Antth    = -sintheta; %dipole radiation pattern

%Projection in the usual rectangular coordinates
Sx8th = (((-sb).^2+(1-(-sb).^2).*ca).*(Antth.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(Antth.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*Antth));
Sy8th = ((-sb.*cb.*(1-ca)).*(Antth.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(Antth.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*Antth));
Sz8th = ((cb.*sa).*(Antth.*costheta.*cosphi)+(sb.*sa).*(Antth.*costheta.*sinphi)+ca.*(-sintheta.*Antth));
AZIMC=acos((POS(:,1)-X_1)./distxy);
AZIMS=asin((POS(:,2)-Y_1)./distxy);
azim=AZIMC.*sign(AZIMS);
elev=asin((POS(:,3)-Z_1)./dist);

