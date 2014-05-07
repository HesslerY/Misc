%CIR.m Computes the channel impulse response at the rectangular coordinates X_1,Y_1,Z_1 for a given sources matrix POS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%           CHANNEL IMPULSE RESPONSE FUNCTION V 2.0             %
%                                                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sx,Sy,Sz,t,azim,elev]=CIRvectca(X_1,Y_1,Z_1)

global Lt c R POS

%Channel Impulse Response computation, we compute 8 CIR corresponding to
%each eighth of the system


%1/4
%CIR 1/8, images in the space defined by x>=0, y>=0, z>=0
[Sx1,Sy1,Sz1,dist1,azim1,elev1] = CIR8thvect(X_1,Y_1,Z_1);

Sx =Sx1;
Sy =Sy1; 
Sz =Sz1;
t =dist1/c;
azim=azim1;
elev=elev1;

