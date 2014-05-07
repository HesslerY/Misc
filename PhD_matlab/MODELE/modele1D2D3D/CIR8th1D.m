
%CIR8th.m computes the CIR of a given 8th of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%        CHANNEL IMPULSE RESPONSE Eighth FUNCTION V 2.0         %
%                                                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sz8th]=CIR8th1D(X_1)

global Lt c R N POS

Sx8th = zeros(1,N);
Sy8th = zeros(1,N);
Sz8th = zeros(1,N);

for j=1:1:length(POS) %loop over the image-currents... a vector version exists and will be released soon
    
    DX = X_1-POS(j,1);
    %DY = Y_1-POS(j,2);
    %DZ = Z_1-POS(j,3);
    
    dist = sqrt(DX^2);%+DY^2+DZ^2);
    zl = round((N-1)*dist/c/Lt);
   
    if zl<N
Vz=1./DX*R^(POS(j,2))*(-1)^POS(j,2);
        
        

        Sz8th(zl+1) = Sz8th(zl+1)+Vz;
        
    end
end