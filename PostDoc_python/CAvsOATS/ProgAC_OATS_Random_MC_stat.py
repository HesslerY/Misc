#!/usr/bin/env python


from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})

c = 299792458
R = 10
f = array(arange(1e6,3e9+1e6,1e6))
M=500
#measurement points
phi=2*pi*rand(M,1)
theta=arccos(2*rand(M,1)-1)

Macth=zeros(len(f))
Macph=zeros(len(f))
Moatsth=zeros(len(f))
Moatsph=zeros(len(f))

#for o in range(0,MC):
#	print o
Ethac=zeros((M,len(f)),'complex')
Ephac=zeros((M,len(f)),'complex')
Erac=zeros((M,len(f)),'complex')
Ethoats=zeros((M,len(f)),'complex')
Ephoats=zeros((M,len(f)),'complex')
Eroats=zeros((M,len(f)),'complex')

#R_eut=.25 #m
#n=2#number of dipoles
#I=zeros((n,7))
#theta_eut=arccos(2*rand(n,1)-1)
#phi_eut=2*pi*rand(n,1)
#x=R_eut*cos(phi_eut)*sin(theta_eut)
#y=R_eut*sin(phi_eut)*sin(theta_eut)
#z=R_eut*cos(theta_eut)
#tilt=arccos(2*rand(n,1)-1)
#azimut=2*pi*rand(n,1)
#ld=.1      
#amplitude=rand(n,1)
#phas=pi*rand(n,1)

R_eut=.5 #m
n=9#number of dipoles
I=zeros((n,7))
theta_eut=arccos(2*rand(n,1)-1)
phi_eut=2*pi*rand(n,1)
x=([[-4./8],[-3./8],[-2./8],[-1.8],[0],[1./8],[2./8],[3./8],[4./8]])#array(linspace(1,n,n)/5,1)#R_eut*cos(phi_eut)*sin(theta_eut)
y=zeros((n,1))#R_eut*sin(phi_eut)*sin(theta_eut)
z=zeros((n,1))#R_eut*cos(theta_eut)
tilt=arccos(2*rand(n,1)-1)
azimut=zeros((n,1))#2*pi*rand(n,1)
ld=.1      
amplitude=ones((n,1))
phas=zeros((n,1))#2*pi*rand(n,1)



I=concatenate((x,y,z,tilt,azimut,amplitude,phas), axis=1)
#Free space caracterisation (perfect anechoic chamber)
for i in range(0,M):
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
	#Erac[i,:] = Exx*sin(TH[i,j])*cos(PH[i,j])+Eyy*sin(theta[i,j])*sin(PH[i,j])+Ezz*cos(TH[i,j])	

#OATS

h=1
I1=concatenate((x,y,z+h,tilt,azimut,amplitude,phas), axis=1)
I2=concatenate((x,y,-(z+h),tilt,azimut+pi,amplitude,phas), axis=1)
Ioats=vstack((I1,I2))
#Free space caracterisation (perfect anechoic chamber)
for i in range(0,M):
	X=R*cos(phi[i])*sin(theta[i])
	Y=R*sin(phi[i])*sin(theta[i])
	Z=R*cos(theta[i])
	DX = X-Ioats[:,0]
	DY = Y-Ioats[:,1]
	DZ = Z-Ioats[:,2]
	dist = sqrt(DX**2+DY**2+DZ**2)
	dp=tile(dist, (len(f),1))
	fp=tile(f,(len(dist),1))
	phaseI=tile(Ioats[:,6],(len(f),1))
	phase=2*pi*dp*fp.T/c+phaseI
	ca    = cos(Ioats[:,3])
	sa    = sin(Ioats[:,3])
	cb    = cos(Ioats[:,4])
	sb    = sin(Ioats[:,4])
	distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
	disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
	distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
	DXY=sqrt(DX**2+DY**2)
	distxy = sqrt(distx**2+disty**2)
	costheta = distz/dist
	sintheta = distxy/dist
	cosphi   = distx/distxy
	sinphi   = disty/distxy
	L =tile(Ioats[:,5],(len(f),1))*1/dp*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
	Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
	Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
	Ezz = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
	Ethoats[i,:]= Exx*cos(theta[i])*cos(phi[i])+Eyy*cos(theta[i])*sin(phi[i])-Ezz*sin(theta[i])
	Ephoats[i,:]= -Exx*sin(phi[i])+Eyy*cos(phi[i])
	#Eroats[j,i,:] = Exx*sin(TH[i,j])*cos(PH[i,j])+Eyy*sin(TH[i,j])*sin(PH[i,j])+Ezz*cos(TH[i,j])


Macth=abs(Ethac).max(axis=0)
Macph=abs(Ephac).max(axis=0)

Moatsth=abs(Ethoats).max(axis=0)
Moatsph=abs(Ephoats).max(axis=0)

macth=abs(Ethac).mean(axis=0)
macph=abs(Ephac).mean(axis=0)

moatsth=abs(Ethoats).mean(axis=0)
moatsph=abs(Ephoats).mean(axis=0)

Pmoats=(abs(Ephoats)**2+abs(Ethoats)**2).mean(axis=0)
Pmac=(abs(Ephac)**2+abs(Ethac)**2).mean(axis=0)

PMoats=(abs(Ephoats)**2+abs(Ethoats)**2).max(axis=0)
PMac=(abs(Ephac)**2+abs(Ethac)**2).max(axis=0)

figure(1)
plot(f,Macth,'r-',f,Moatsth,'b-',f,Macph,'r--',f,Moatsph,'b--')
legend((r'$\max(E_\theta)$ (AC)',r'$\max(E_\theta)$ (OATS)',r'$\max(E_{\phi})$ (AC)',r'$\max(E_{\phi})$ (OATS)'),loc=2)
grid('on')

figure(2)
plot(f,macth,'r-',f,moatsth,'b-',f,macph,'r--',f,moatsph,'b--')
legend((r'$<E_\theta>$ (AC)',r'$<E_\theta>$ (OATS)',r'$<E_{\phi}>$ (AC)',r'$<E_{\phi}>$ (OATS)'),loc=2)
grid('on')


figure(3)
plot(f,Moatsth/Macth,f,Moatsph/Macph)
legend((r'$\max(E_{\theta_{oats}})/\max(E_{\theta_{ac}})$',r'$\max(E_{\phi_{oats}})/\max(E_{\phi_{ac}})$'),loc=2)
grid('on')

figure(4)
plot(f,moatsth/macth,f,moatsph/macph)
legend((r'$<E_{\theta_{oats}}>/<E_{\theta_{ac}}>$',r'$<E_{\phi_{oats}}>/<E_{\phi_{ac}}>$'),loc=2)
grid('on')

figure(5)
plot(f,Macth/macth,f,Macph/macph)
legend((r'$\max(E_{\theta_{ac}})/<E_{\theta_{ac}}>$',r'$\max(E_{\phi_{ac}})/<E_{\phi_{ac}}>$'),loc=2)
grid('on')

figure(6)
plot(f,Moatsth/moatsth,f,Moatsph/moatsph)
legend((r'$\max(E_{\theta_{oats}})/<E_{\theta_{oats}}>$',r'$\max(E_{\phi_{oats}})/<E_{\phi_{oats}}>$'),loc=2)
grid('on')

figure(7)

plot(f,Pmoats/Pmac,f,PMoats/PMac)


show()
