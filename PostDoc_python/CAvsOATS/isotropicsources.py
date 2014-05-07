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
f = array(arange(100e6,10e9+10e6,100e6))
#f = array([1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])*1e6
np=360
nt=180
dtheta = pi/nt
dphi = (2*pi)/np
#measurement points
phi=linspace(0,2*pi,np)
theta=linspace(0,pi,nt)#arccos(2*rand(M,1)-1)
Dac=zeros(len(f))
Dacp=zeros(len(f))
#MC=500
#Dac=zeros((MC,len(f)))
#Doats=zeros((MC,len(f)))
start = time.time()
#for o in range(0,MC):  
E=zeros((len(phi),len(theta),len(f)))

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
		L =I[:,5]*1/dp
		E[i,j,:]= abs((L*exp(1J*phase)).sum(axis=1))

for u in range(0,len(f)):
	Fa=E[:,:,u]**2
	Dacp[u]=Fa.max()/Fa.mean()

ka=2*pi*f/c*R_eut
k=array(arange(1,110,1))
Dmaxth=(0.577+log(4*k**2+8*k)+1/(8*k**2+16*k))

figure(1)
plot(ka,Dacp,k,Dmaxth)


show()

