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
f = array(arange(5e6,3e9+10e6,5e6))
#f = array([1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])*1e6

#measurement points
M=100;

MC=30
Pmac=zeros((MC,len(f)))
Pmoats=zeros((MC,len(f)))
PMac=zeros((MC,len(f)))
PMoats=zeros((MC,len(f)))
start = time.time()
for o in range(0,MC):	
	PH=2*pi*rand(M,1)
	TH=arccos(2*rand(M,1)-1)	
	Ethac=zeros((M,len(f)),'complex')
	Ephac=zeros((M,len(f)),'complex')
	#Erac=zeros((len(theta),len(phi),len(f)),'complex')
	Ethoats=zeros((M,len(f)),'complex')
	Ephoats=zeros((M,len(f)),'complex')
	#Eroats=zeros((len(theta),len(phi),len(f)),'complex')
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
		X=R*cos(PH[i])*sin(TH[i])
		Y=R*sin(PH[i])*sin(TH[i])
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
	I1=concatenate((x,y,z+h,tilt,azimut,amplitude,phas), axis=1)
	I2=concatenate((x,y,-(z+h),tilt,azimut+pi,amplitude,phas), axis=1)
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
		#Eroats[i,:] = Exx*sin(TH[i])*cos(PH[i])+Eyy*sin(TH[i])*sin(PH[i])+Ezz*cos(TH[i])
	for u in range(0,len(f)):
		Pac=abs(Ephac[:,u])**2+abs(Ethac[:,u])**2
		Poats=abs(Ephoats[:,u])**2+abs(Ethoats[:,u])**2
		Pmac[o,u]=mean(Pac)
		Pmoats[o,u]=mean(Poats)
		PMac[o,u]=Pac.max()/mean(Pac)
		PMoats[o,u]=Poats.max()/mean(Poats)
	elapsed = (time.time() - start)
	print ('%d, t= %2.2f s' %(o,elapsed))
	

figure(1)
subplot(211)
plot(f,mean(Pmac,axis=0),f,mean(Pmoats,axis=0))
subplot(212)
plot(f,mean(PMac,axis=0),f,mean(PMoats,axis=0))

show()









#savetxt('Doats_500_n30_R25cm.txt', Doats)
#savetxt('Dac_500_n30_R25cm.txt', Dac)
#savetxt('f_500_n30_R25cm.txt', f)

#ka=2*pi*f/c*R_eut
#k=array(arange(1,200,1))
#Dmaxth=1./2*(0.577+log(4*k**2+8*k)+1/(8*k**2+16*k))
#
#figure(1)
#plot(f/1e6,mean(Dac,axis=0),f/1e6,mean(Doats,axis=0))#,f/1e6,Dmaxth)
#grid('on')
#xlim((0, 3000))
#xlabel(r'$f$ [MHz]')
#ylabel(r'$<D>$')
#legend((r'$<D>_{\textrm{far}}$',r'$<D>_{\textrm{oats}}$'),loc=2)
#
#figure(2)
#plot(f/1e6,mean(Doats,axis=0)/mean(Dac,axis=0))#,f/1e6,Dmaxth)
#grid('on')
#xlim((0, 3000))
#xlabel(r'$f$ [MHz]')
#ylabel(r'$<D>$')
##legend((r'$\frac{<D>_{\textrm{oats}}}{<D>_{\textrm{oats}}}$'),loc=1)
#title(r'$<D>_{\textrm{oats}}/<D>_{\textrm{ac}}$')
#
#figure(4)
#plot(ka*4,mean(Dac,axis=0),'b--',ka*4,mean(Doats,axis=0),'r--',k,Dmaxth,'k')
#grid('on')
#xlim((0, max(ka*4)))
#legend((r'$<D>_{\textrm{far}}$ with $a=4R_\textrm{eut}$',r'$<D>_{\textrm{oats}}$ with $a=4R_\textrm{eut}$',r'$<D_{\max}>_\textrm{th}$'),loc=2)
#xlabel(r'$ka$')
#ylabel(r'$<D>$')
#
#figure(5)
#plot(ka*4,mean(Dac,axis=0),ka*4*1.5,mean(Doats,axis=0),k,Dmaxth)
#grid('on')
#xlim((0, max(ka*6)))
#show()
#

