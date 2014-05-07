from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams
import time
import os
import sys
files = []

c = 299792458
R = 10.
f = array(arange(100e6,10e9+100e6,100e6))
#f = array([1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])*1e6
np=360
nt=180
dtheta = pi/nt
dphi = (2*pi)/np
#measurement points
phi=linspace(0,2*pi,np)
theta=linspace(0,pi,nt)#arccos(2*rand(M,1)-1)
TH,PH=meshgrid(theta,phi)
Dac=zeros(len(f))
Dacp=zeros(len(f))
#MC=500
#Dac=zeros((MC,len(f)))
#Doats=zeros((MC,len(f)))
start = time.time()
#for o in range(0,MC):  
Ethac=zeros((len(phi),len(theta),len(f)),'complex')
Ephac=zeros((len(phi),len(theta),len(f)),'complex')
P=zeros((len(phi),len(theta),len(f)))
#Eroats=zeros((len(z),len(phi),len(f)),'complex')
#FAR
R_eut=.5 #m
n=30#number of dipoles
I=zeros((n,7))
theta_eut=arccos(2*rand(n,1)-1)
phi_eut=2*pi*rand(n,1)
xx=R_eut*cos(phi_eut)*sin(theta_eut)
yy=R_eut*sin(phi_eut)*sin(theta_eut)
zz=R_eut*cos(theta_eut)
tilt=arccos(2*rand(n,1)-1)
azimut=2*pi*rand(n,1)
amplitude=rand(n,1)
phas=2*pi*rand(n,1)
I=concatenate((xx,yy,zz,tilt,azimut,amplitude,phas), axis=1)

#isotropic
for i in range(0,len(phi)):
	for j in range(0,len(theta)):
		X=R*cos(phi[i])*sin(theta[j])
		Y=R*sin(phi[i])*sin(theta[j])
		Z=R*cos(theta[j])
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
                L =tile(I[:,5],(len(f),1))*1/dp#*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
                Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
                Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
                Ezz = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
                Ethac[i,j,:]= Exx*cos(theta[j])*cos(phi[i])+Eyy*cos(theta[j])*sin(phi[i])-Ezz*sin(theta[j])
                Ephac[i,j,:]= -Exx*sin(phi[i])+Eyy*cos(phi[i])
                P[i,j,:]=abs(Exx*cos(theta[j])*cos(phi[i])+Eyy*cos(theta[j])*sin(phi[i])-Ezz*sin(theta[j]))**2+abs(-Exx*sin(phi[i])+Eyy*cos(phi[i]))**2

for u in range(0,len(f)):
	Fa2ac=abs(Ephac[:,:,u])**2+abs(Ethac[:,:,u])**2
	Faac=Fa2ac/Fa2ac.max()
	omegaac = (Faac*sin(TH)*dtheta*dphi).sum()
	Dac[u] = 4*pi/omegaac
	Dacp[u]=P[:,:,u].max()/P[:,:,u].mean()

ka=2*pi*f/c*2*R_eut
k=array(arange(1,200,1))
Dmaxth=1./2*(0.577+log(4*k**2+8*k)+1/(8*k**2+16*k))

figure(1)
plot(ka,Dacp,ka,Dac,k,Dmaxth)


show()

