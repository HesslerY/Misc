#!/usr/bin/env python


from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams
import time

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})

c = 299792458
R = 10
f = array(arange(5e6,2e9+5e6,5e6))
M=100;
#measurement points




MC=30
Powac=zeros((MC,len(f)))
Powoats=zeros((MC,len(f)))
start = time.time()
for o in range(0,MC):	
	PH=2*pi*rand(M,1)
	TH=arccos(2*rand(M,1)-1)
	Ethac=zeros((len(PH),len(f)),'complex')
	Ephac=zeros((len(PH),len(f)),'complex')
	Ethoats=zeros((len(PH),len(f)),'complex')
	Ephoats=zeros((len(PH),len(f)),'complex')
	R_eut=.25 #m
	n=20#number of dipoles
	I=zeros((n,7))
	theta_eut=arccos(2*rand(n,1)-1)
	phi_eut=2*pi*rand(n,1)
	x=R_eut*cos(phi_eut)*sin(theta_eut)
	y=R_eut*sin(phi_eut)*sin(theta_eut)
	z=R_eut*cos(theta_eut)
	tilt=arccos(2*rand(n,1)-1)
	azimut=2*pi*rand(n,1)
	ld=.1      
	amplitude=rand(n,1)
	phas=2*pi*rand(n,1)
	I=concatenate((x,y,z,tilt,azimut,amplitude,phas), axis=1)
	#Free space caracterisation (perfect anechoic chamber)
	for i in range(0,len(PH)):
		X=R*cos(PH[i])*sin([i])
		Y=R*sin(PH[i])*sin([i])
		Z=R*cos(TH[i])
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
		Ethac[i,:]= Exx*cos(TH[i])*cos(PH[i])+Eyy*cos(TH[i])*sin(PH[i])-Ezz*sin(TH[i])
		Ephac[i,:]= -Exx*sin(PH[i])+Eyy*cos(PH[i])
		#Erac[i,:] = Exx*sin(TH[i])*cos(PH[i])+Eyy*sin(TH[i])*sin(PH[i])+Ezz*cos(TH[i])	
	#OATS
	h=1
	I1=concatenate((x,y,z,tilt,azimut,amplitude,phas), axis=1)
	I2=concatenate((x,y,-(z+2*h),tilt,azimut+pi,amplitude,phas), axis=1)
	Ioats=vstack((I1,I2))
	#Free space caracterisation (perfect anechoic chamber)
	for i in range(0,len(PH)):
		X=R*cos(PH[i])*sin(TH[i])
		Y=R*sin(PH[i])*sin(TH[i])
		Z=R*cos(TH[i])
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
		Ethoats[i,:]= Exx*cos(TH[i])*cos(PH[i])+Eyy*cos(TH[i])*sin(PH[i])-Ezz*sin(TH[i])
		Ephoats[i,:]= -Exx*sin(PH[i])+Eyy*cos(PH[i])
		#Eroats[i,:] = Exx*sin(TH[i])*cos(PH[i])+Eyy*sin(TH[i][i])*sin(PH[i])+Ezz*cos(TH[i])
	for u in range(0,len(f)):
		Pac=abs(Ethac[:,u])**2+abs(Ephac[:,u])**2
		Poats=abs(Ethoats[:,u])**2+abs(Ephoats[:,u])**2
		Powac[o,u] = Pac.max()/mean(Pac)
		Powoats[o,u] = Poats.max()/mean(Poats)
	elapsed = (time.time() - start)
	print ('%d, t= %2.2f s' %(o,elapsed))
	
#savetxt('D1Doats_100_n30_R25cm.txt', Doats)
#savetxt('D1Dac_100_n30_R25cm.txt', Dac)
#savetxt('f1D_100_n30_R25cm.txt', f)
#

figure(1)
plot(f/1e6,mean(Powac,axis=0),f/1e6,mean(Powoats,axis=0))#,f/1e6,Dmaxth)
grid('on')
xlim((0, 3000))
xlabel(r'$f$ [MHz]')
ylabel(r'$<D>$')
legend((r'$<D>_{\textrm{far}}$',r'$<D>_{\textrm{oats}}$'),loc=2)

figure(2)
plot(f/1e6,mean(Powoats,axis=0)/mean(Powac,axis=0))#,f/1e6,Dmaxth)
grid('on')
xlim((0, 3000))
xlabel(r'$f$ [MHz]')
ylabel(r'$<D>$')
#legend((r'$\frac{<D>_{\textrm{oats}}}{<D>_{\textrm{oats}}}$'),loc=1)
title(r'$<D>_{\textrm{oats}}/<D>_{\textrm{ac}}$')

show()
