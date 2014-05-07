from numpy import *
from numpy.random import *
from pylab import *

c = 299792458
R = 10
f = array(arange(50e6,1e9,50e6))

#measurement points
M = 1000
phi=2*pi*rand(M,1)
theta=arccos(2*rand(M,1)-1)

#EUT
R_eut=.3 #m
n=20 #number of dipoles
I=zeros((n,7))
theta_eut=arccos(2*rand(n,1)-1)
phi_eut=2*pi*rand(n,1)
x=R_eut*cos(phi_eut)*sin(theta_eut)
y=R_eut*sin(phi_eut)*sin(theta_eut)
z=R_eut*cos(theta_eut)
tilt=arccos(2*rand(n,1)-1)
azimut=2*pi*rand(n,1)
ld=.1      
amplitude=ones((n,1))
phas=2*pi*rand(n,1)
I=concatenate((x,y,z,tilt,azimut,amplitude,phas), axis=1)

Ethac=zeros((len(theta),len(f)),'complex')
Ephac=zeros((len(theta),len(f)),'complex')
Erac=zeros((len(theta),len(f)),'complex')
#Free space caracterisation (perfect anechoic chamber)
for i in range(0,len(theta)):
	X=R*cos(phi[i])*sin(theta[i])
	Y=R*sin(phi[i])*sin(theta[i])
	Z=R*cos(theta[i])
	DX = X-I[:,0]
	DY = Y-I[:,1]
	DZ = Z-I[:,2]
	dist = sqrt(DX**2+DY**2+DZ**2)
	dp=tile(dist, (len(f),1))
	fp=tile(f,(len(dist),1))
	phaseI=tile(I[:,6],(len(f),1))
	phase=2*pi*dp*fp.T/c+phaseI
	ca    = cos(I[:,3])
	sa    = sin(I[:,3])
	cb    = cos(I[:,4])
	sb    = sin(I[:,4])
	distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
	disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
	distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
	DXY=sqrt(DX**2+DY**2)
	distxy = sqrt(distx**2+disty**2)
	costheta = distz/dist
	sintheta = distxy/dist
	cosphi   = distx/distxy
	sinphi   = disty/distxy
	L =tile(I[:,5],(len(f),1))*1/dp*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
	Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
	Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
	Ezz = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
	Ethac[i,:]= Exx*cos(theta[i])*cos(phi[i])+Eyy*cos(theta[i])*sin(phi[i])-Ezz*sin(theta[i])
	Ephac[i,:]= -Exx*sin(phi[i])+Eyy*cos(phi[i])
	Erac[i,:] = Exx*sin(theta[i])*cos(phi[i])+Eyy*sin(theta[i])*sin(phi[i])+Ezz*cos(theta[i])

for i in range(1,len(f)):
	figure(i)
	subplot(131)
	hist(abs(Ethac[:,i]),30)
	subplot(132)
	hist(abs(Ephac[:,i]),30)
	subplot(133)
	hist(abs(Erac[:,i]),30)

show()