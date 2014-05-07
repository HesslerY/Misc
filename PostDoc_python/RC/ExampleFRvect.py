#!/usr/bin/python

"""Frequency response in an EM reverberation chamber
Visit http://imagetheorymodel.manuamador.fr/ for the
latest version and more information
"""

__author__ = "Emmanuel Amador (emmanuel.amador@insa-rennes.fr)"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2011/06/12$"
__copyright__ = "Copyright (c) 2011 E. Amador"
__license__ = "NDA"



from numpy import *
from pylab import *

def IC(Tm,l,p,h,X,Y,Z,tilt,azimut): #Computes images' postions and image's orientation
	dims = array([l,p,h])
	c = 299792458
	Tp = Tm+3./c*dims.max()
	dmax=c*Tp
	order = round(dmax/dims.min())
	Memo = pi/6*(Tp*c)**3/l/p/h*6*8
	POSp = array([X,Y,Z, 0, tilt, azimut]);
	#1D
	for i in range(1,int(order)):
		POSp = vstack([POSp,[2*i*l-X, Y, Z, abs(2*i)-1, pi-tilt, 2*pi-azimut],[2*i*l+X, Y, Z, abs(2*i), tilt, azimut]])
	POSi=POSp.copy();
	POSi[:,1] = array([p*ones((1,len(POSp)))-POSi[:,1]])
	POSi[:,4] = array([pi*ones((1,len(POSp)))-POSi[:,4]])
	POSi[:,5] = array([pi*ones((1,len(POSp)))-POSi[:,5]])
	
	#2D
	POSpp = POSp.copy();
	POSpp[:,1] = array([p*ones((1,len(POSp)))]+POSp[:,1])
	POSpp[:,3] = array([ones((1,len(POSp)))]+POSp[:,3])
	POSii=POSi.copy()
	POSii[:,1] = array([p*ones((1,len(POSi)))]+POSi[:,1])
	POSii[:,3] = array([ones((1,len(POSi)))]+POSi[:,3])
	POS1 = POSpp.copy()
	POS1 = concatenate((POS1,POSii),axis=0)
	dist = (POS1[:,0]**2+POS1[:,1]**2+POS1[:,2]**2)**.5
	POS1 = concatenate((POS1,dist.reshape(-1,1)),axis=1)
	POS1 = POS1[POS1[:,6].argsort(),]
	if POS1[:,6].max()>c*Tp:
		U = where(POS1[:,6]>c*Tp)
		POS1=delete(POS1, s_[min(min(U)):], axis=0)
	
	POS1 = delete(POS1, s_[6], axis=1)
	POSP=POS1.copy()
	
	for jj in range (2,int(order),2):
		POSpp = POSp.copy()
		POSpp[:,1] = array([jj*p*ones((1,len(POSp)))]+POSp[:,1])
		POSpp[:,3] = array([abs(jj)*ones((1,len(POSp)))]+POSp[:,3])
		POSii=POSi.copy()
		POSii[:,1] = array([(jj+1)*p*ones((1,len(POSi)))]+POSi[:,1])
		POSii[:,3] = array([abs(jj+1)*ones((1,len(POSi)))]+POSi[:,3])
		POS2 = POSpp.copy()
		POS2 = concatenate((POS2,POSii),axis=0)
		dist = (POS2[:,0]**2+POS2[:,1]**2+POS2[:,2]**2)**.5
		POS2 = concatenate((POS2,dist.reshape(-1,1)),axis=1)
		POS2 = POS2[POS2[:,6].argsort(),]
		if POS2[:,6].max()>c*Tp:
			V = where(POS2[:,6]>c*Tp)
			POS2=delete(POS2, s_[min(min(V)):], axis=0)
		POS2 = delete(POS2, s_[6], axis=1)
		POSP = concatenate((POSP,POS2),axis=0)
		
	
	POSI = POSP.copy();
	POSI[:,2] = array([h*ones((1,len(POSI)))-POSI[:,2]])
	POSI[:,5] = mod(array([pi*ones((1,len(POSI)))+POSI[:,5]]),2*pi)
	
	#3D	
	POSPP = POSP.copy()
	POSII=POSI.copy()
	POSII[:,2] = array([h*ones((1,len(POSI)))]+POSI[:,2])
	POSII[:,3] = array([ones((1,len(POSI)))]+POSI[:,3])
	POS3 = POSPP.copy()
	POS3 = concatenate((POS3,POSII),axis=0)
	dist = (POS3[:,0]**2+POS3[:,1]**2+POS3[:,2]**2)**.5
	POS3 = concatenate((POS3,dist.reshape(-1,1)),axis=1)
	POS3 = POS3[POS3[:,6].argsort(),]
	if POS3[:,6].max()>c*Tp:
		W = where(POS3[:,6]>c*Tp)
		POS3=delete(POS3, s_[min(min(W)):], axis=0)	
	POS3 = delete(POS3, s_[6], axis=1)
	POS=POS3.copy()		
	
	for k in range (2,int(order),2):
		POSPP = POSP.copy()
		POSPP[:,2] = array([k*h*ones((1,len(POSP)))]+POSP[:,2])
		POSPP[:,3] = array([(k)*ones((1,len(POSP)))]+POSP[:,3])
		POSII=POSI.copy()
		POSII[:,2] = array([(k+1)*h*ones((1,len(POSI)))]+POSI[:,2])
		POSII[:,3] = array([(k+1)*ones((1,len(POSI)))]+POSI[:,3])
		POS4 = POSPP.copy()
		POS4 = concatenate((POS4,POSII),axis=0)
		dist = (POS4[:,0]**2+POS4[:,1]**2+POS4[:,2]**2)**.5
		POS4 = concatenate((POS4,dist.reshape(-1,1)),axis=1)
		POS4 = POS4[POS4[:,6].argsort(),]
		if POS4[:,6].max()>c*Tp:
			Q = where(POS4[:,6]>c*Tp)
			POS4=delete(POS4, s_[min(min(Q)):], axis=0)
		POS4 = delete(POS4, s_[6], axis=1)
		POS = concatenate((POS,POS4),axis=0)
	
	POS=vstack((POSp[0,:],POS))
	return POS



def FR(a,b,c):
	#1/8
	Fx1,Fy1,Fz1 = FR8th(a,b,c)
	2/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx2,Fy2,Fz2 = FR8th(a,b,c)
	POS[:,2]=-POS[:,2]
 	3/8
	POS[:,1]=-POS[:,1]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	Fx3,Fy3,Fz3 = FR8th(a,b,c)
	4/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx4,Fy4,Fz4 = FR8th(a,b,c)
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	5/8
	POS[:,3]=array([POS[:,3]-ones((1,len(POS)))])
	POS[:,0]=-POS[:,0]
	Fx5,Fy5,Fz5 = FR8th(a,b,c)
	6/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx6,Fy6,Fz6 = FR8th(a,b,c)
	POS[:,2]=-POS[:,2]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	7/8
	POS[:,1]=-POS[:,1]
	POS[:,5]=-POS[:,5]
	Fx7,Fy7,Fz7 = FR8th(a,b,c)
	8/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Fx8,Fy8,Fz8 = FR8th(a,b,c)
	POS[:,0]=-POS[:,0]
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	POS[:,3]=array([-3.*ones((1,len(POS)))+POS[:,3]])
	POS[:,5]=mod(POS[:,5],2*pi)
	Fx=array([Fx1]+[Fx2]+[Fx3]+[Fx4]+[Fx5]+[Fx6]+[Fx7]+[Fx8]).sum(axis=0)
	Fy=array([Fy1]+[Fy2]+[Fy3]+[Fy4]+[Fy5]+[Fy6]+[Fy7]+[Fy8]).sum(axis=0)
	Fz=array([Fz1]+[Fz2]+[Fz3]+[Fz4]+[Fz5]+[Fz6]+[Fz7]+[Fz8]).sum(axis=0)
	return Fx,Fy,Fz


def FR8th(r,s,t):
	Ex8th=zeros((len(r),len(f)),'complex')
	Ey8th=zeros((len(r),len(f)),'complex')
	Ez8th=zeros((len(r),len(f)),'complex')
	for i in range(0,len(r)):
		DX = r[i]-POS[:,0]
		DY = s[i]-POS[:,1]
		DZ = t[i]-POS[:,2]
		dist = sqrt(DX**2+DY**2+DZ**2)
		dp=tile(dist, (len(f),1))
		fp=tile(f,(len(dist),1))
		phase=2*pi*dp*fp.T/c
		ca    = cos(POS[:,4])
		sa    = sin(POS[:,4])
		cb    = cos(POS[:,5])
		sb    = sin(POS[:,5])
		distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
		disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
		distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
		DXY=sqrt(DX**2+DY**2)
		distxy = sqrt(distx**2+disty**2)
		costheta = distz/dist
		sintheta = distxy/dist
		cosphi   = distx/distxy
		sinphi   = disty/distxy
		L =tile(R**(POS[:,3]),(len(f),1))*1/dp #Amplitude & free space attenuation
		Ex8th[i,:] = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Ey8th[i,:] = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Ez8th[i,:] = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
	return Ex8th,Ey8th,Ez8th

c = 299792458

#Image Creation
Lt=1e-6 #length of the time window in s

#dimensions of the reverb chamber in m
l=8.7
p=3.7
h=2.9

#position of the emitter and angular orientation		
X = 1.
Y = 2.
Z = 1.
tilt = pi/2-math.acos(sqrt(2./3));
azimut = pi/4;

POS=IC(Lt,l,p,h,X,Y,Z,tilt,azimut)

#Frequency response calculation
f = array(arange(10e6,150e6,.5e6))

R=0.998 #loss coefficient

#measurement point
X_1=array([4.5])
Y_1=array([3])
Z_1=array([1])

Fx,Fy,Fz=FR(X_1,Y_1,Z_1)

for i in range(0,len(Fz[:,0])):
	figure(i)
	subplot(311)
	l = plot(f/1e6,20*log10(abs(Fx[i,:])))
	grid(True)
	title('Frequency response')
	ylabel('$E_x$')	
	subplot(312)
	l = plot(f/1e6,20*log10(abs(Fy[i,:])))
	grid(True)
	ylabel('$E_y$')
	subplot(313)
	l = plot(f/1e6,20*log10(abs(Fz[i,:])))
	ylabel('$E_z$')
	grid(True)
	xlabel('frequency (MHz)')

show()
