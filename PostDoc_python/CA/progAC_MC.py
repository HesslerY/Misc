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
f = array(arange(1e6,1e9+1e6,1e6))

#measurement points
N=20
M = 50
MC=100
phi=2*pi*rand(M,1)
theta=arccos(2*rand(M,1)-1)
D=zeros((N,MC,len(f)))
Ethac=zeros((20,len(theta),len(f)),'complex')
Ephac=zeros((20,len(theta),len(f)),'complex')
Erac=zeros((20,len(theta),len(f)),'complex')

for o in range(0,MC,1):
	#EUT
	R_eut=.3 #m
	for n in range(1,N+1,1): #number of dipoles
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
			Ethac[n-1,i,:]= Exx*cos(theta[i])*cos(phi[i])+Eyy*cos(theta[i])*sin(phi[i])-Ezz*sin(theta[i])
			Ephac[n-1,i,:]= -Exx*sin(phi[i])+Eyy*cos(phi[i])
			Erac[n-1,i,:] = Exx*sin(theta[i])*cos(phi[i])+Eyy*sin(theta[i])*sin(phi[i])+Ezz*cos(theta[i])
	
	
	for i in range(0,N,1):
		A=vstack((Ethac[i,:,:],Ephac[i,:,:]))
		D[i,o,:]=mean(abs(A),axis=0)

moy=zeros((n,len(f)))
per2=zeros((n,len(f)))
per9=zeros((n,len(f)))
for i in range(0,n,1):
	for j in range(0,len(f),1):
		moy[i,j]=mean(D[i,:,j],axis=0)
		per2[i,j]=prctile(D[i,:,j], p=(2.5))
		per9[i,j]=prctile(D[i,:,j], p=(97.5))

nn=linspace(1,20,20);
Erth=377*2*pi*1e9/c*ld/16/R*sqrt(nn*pi/2)/2;


plot(nn,Erth*pi/3,'r--',nn,moy[:,19],'k-',nn,per2[:,19],'k--',nn,per9[:,19],'k--')
grid('on')
legend((r'$<E_R(N)> \textrm{(6)}$',r'$<E_R(N)>$', r'$\textrm{95\% CI of }<E_R(N)>$'),loc=2 )
ylabel(r'$\textrm{(V/m)}$')
xlabel('$N$')
show()


