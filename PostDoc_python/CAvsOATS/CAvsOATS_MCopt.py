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
f = array(arange(1e6,5e9+50e6,50e6))
#f = array([1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])*1e6
np=90
nt=45
dtheta = pi/nt
dphi = (2*pi)/np
#measurement points
phi=linspace(0,2*pi,np)
theta=linspace(0,pi,nt)#arccos(2*rand(M,1)-1)
TH,PH=meshgrid(theta,phi)

#MC=500
#Dac=zeros((MC,len(f)))
#Doats=zeros((MC,len(f)))
#for o in range(0,MC):  

#Eroats=zeros((len(z),len(phi),len(f)),'complex')
#FAR

MC=1000
Dac=zeros((MC,len(f)))
Dacp=zeros((MC,len(f)))
Doats=zeros((MC,len(f)))
Doatsp=zeros((MC,len(f)))



for o in range (0,MC):#isotropic
	start = time.time()
	Ethac=zeros((len(phi),len(theta),len(f)),'complex')
	Ephac=zeros((len(phi),len(theta),len(f)),'complex')
	Pac=zeros((len(phi),len(theta),len(f)))
	Ethoats=zeros((len(phi),len(theta),len(f)),'complex')
	Ephoats=zeros((len(phi),len(theta),len(f)),'complex')
	Poats=zeros((len(phi),len(theta),len(f)))
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
	h=1
	I1=concatenate((xx,yy,zz,tilt,azimut,amplitude,phas), axis=1)
	I2=concatenate((xx,yy,-(zz+2*h),tilt,azimut+pi,amplitude,phas), axis=1)
	for i in range(0,len(phi)):
		for j in range(0,len(theta)):
			X=R*cos(phi[i])*sin(theta[j])
			Y=R*sin(phi[i])*sin(theta[j])
			Z=R*cos(theta[j])
			DX = X-I1[:,0]
			DY = Y-I1[:,1]
			DZ = Z-I1[:,2]
			dist = sqrt(DX**2+DY**2+DZ**2)
			dp=tile(dist, (len(f),1))
			fp=tile(f,(len(dist),1))
			phaseI=tile(I1[:,6],(len(f),1))
			phase=2*pi*dp*fp.T/c+phaseI
			ca    = cos(I1[:,3])
			sa    = sin(I1[:,3])
			cb    = cos(I1[:,4])
			sb    = sin(I1[:,4])
			distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
			disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
			distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
			DXY=sqrt(DX**2+DY**2)
			distxy = sqrt(distx**2+disty**2)
			costheta = distz/dist
			sintheta = distxy/dist
			cosphi   = distx/distxy
			sinphi   = disty/distxy
			L =tile(I1[:,5],(len(f),1))*1/dp#*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
			Exx1 = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
			Eyy1 = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
			Ezz1 = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
			Ethac[i,j,:]= Exx1*cos(theta[j])*cos(phi[i])+Eyy1*cos(theta[j])*sin(phi[i])-Ezz1*sin(theta[j])
			Ephac[i,j,:]= -Exx1*sin(phi[i])+Eyy1*cos(phi[i])
			Pac[i,j,:]=abs(Exx1*cos(theta[j])*cos(phi[i])+Eyy1*cos(theta[j])*sin(phi[i])-Ezz1*sin(theta[j]))**2+abs(-Exx1*sin(phi[i])+Eyy1*cos(phi[i]))**2
			DX1 = X-I2[:,0]
			DY1 = Y-I2[:,1]
			DZ1 = Z-I2[:,2]
			dist1 = sqrt(DX1**2+DY1**2+DZ1**2)
			dp1=tile(dist1, (len(f),1))
			fp1=tile(f,(len(dist1),1))
			phaseI1=tile(I2[:,6],(len(f),1))
			phase1=2*pi*dp1*fp1.T/c+phaseI1
			ca1    = cos(I2[:,3])
			sa1    = sin(I2[:,3])
			cb1    = cos(I2[:,4])
			sb1    = sin(I2[:,4])
			distx1 = ((-sb1)**2+(1-(-sb1)**2)*ca1)*DX1+(-sb1*cb1*(1-ca1))*DY1+(cb1*sa1)*DZ1
			disty1 = (-sb1*cb1*(1-ca1))*DX1+((cb1)**2+(1-cb1**2)*ca1)*DY1+(sb1*sa1)*DZ1
			distz1 = (-cb1*sa1)*DX1+(-sb1*sa1)*DY1+ca1*DZ1
			DXY1=sqrt(DX1**2+DY1**2)
			distxy1 = sqrt(distx1**2+disty1**2)
			costheta1 = distz1/dist1
			sintheta1 = distxy1/dist1
			cosphi1   = distx1/distxy1
			sinphi1   = disty1/distxy1
			L1 =tile(I2[:,5],(len(f),1))*1/dp1#*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
			Exx2 =sum(exp(1j*phase1)*L1*tile(((((-sb1)**2+(1-(-sb1)**2)*ca1)*(-sintheta1*costheta1*cosphi1)+(-sb1*cb1*(1-ca1))*(-sintheta1*costheta1*sinphi1)+(-cb1*sa1)*(-sintheta1*(-sintheta1)))),(len(f),1)),axis=1)
			Eyy2 =sum(exp(1j*phase1)*L1*tile((((-sb1*cb1*(1-ca1))*(-sintheta1*costheta1*cosphi1)+((cb1)**2+(1-(cb1)**2)*ca1)*(-sintheta1*costheta1*sinphi1)+(-sb1*sa1)*(-sintheta1*(-sintheta1)))),(len(f),1)),axis=1)
			Ezz2 =sum(exp(1j*phase1)*L1*tile((((cb1*sa1)*(-sintheta1*costheta1*cosphi1)+(sb1*sa1)*(-sintheta1*costheta1*sinphi1)+ca1*(-sintheta1*(-sintheta1)))),(len(f),1)),axis=1)
			Ethoats[i,j,:]= (Exx1.copy()+Exx2)*cos(theta[j])*cos(phi[i])+(Eyy1.copy()+Eyy2)*cos(theta[j])*sin(phi[i])-(Ezz1.copy()+Ezz2)*sin(theta[j])
			Ephoats[i,j,:]= -(Exx1.copy()+Exx2)*sin(phi[i])+(Eyy1.copy()+Eyy2)*cos(phi[i])
			Poats[i,j,:]=abs((Exx1.copy()+Exx2)*cos(theta[j])*cos(phi[i])+(Eyy1.copy()+Eyy2)*cos(theta[j])*sin(phi[i])-(Ezz1.copy()+Ezz2)*sin(theta[j]))**2+abs(-(Exx1.copy()+Exx2)*sin(phi[i])+(Eyy1.copy()+Eyy2)*cos(phi[i]))**2
	for u in range(0,len(f)):
		Fa2ac=abs(Ephac[:,:,u])**2+abs(Ethac[:,:,u])**2
		Faac=Fa2ac/Fa2ac.max()
		omegaac = (Faac*sin(TH)*dtheta*dphi).sum()
		Dac[o,u] = 4*pi/omegaac
		Dacp[o,u]=Pac[:,:,u].max()/Pac[:,:,u].mean()
		Fa2oats=abs(Ephoats[:,:,u])**2+abs(Ethoats[:,:,u])**2
		Faoats=Fa2oats/Fa2oats.max()
		omegaoats = (Faoats*sin(TH)*dtheta*dphi).sum()
		Doats[o,u] = 4*pi/omegaoats
		Doatsp[o,u]=Poats[:,:,u].max()/Poats[:,:,u].mean()
	elapsed = (time.time() - start)
	minutes, secondes = divmod(elapsed*(MC-o-1), 60)
	heures, minutes = divmod(minutes, 60)
	print ('Exp. # %d/%d, dt = %2.2f s, ETA = %d:%02d:%02d' %(o+1,MC,elapsed,heures, minutes, secondes))

savetxt('Doatsp_1000b_n30_R50cm_h1m.txt', Doatsp)
savetxt('Dacp_1000b_n30_R50cm_h1m.txt', Dacp)
savetxt('Doats_1000b_n30_R50cm_h1m.txt', Doats)
savetxt('Dac_1000b_n30_R50cm_h1m.txt', Dac)
savetxt('f_1000b_n30_R50cm_h1m.txt', f)

ka=2*pi*f/c*2*R_eut
k=array(arange(1,110,1))
Dmaxth=1./2*(0.577+log(4*k**2+8*k)+1/(8*k**2+16*k))





figure(1)
#plot(k,Dmaxth,ka,mean(Dacp,axis=0),ka,mean(Dac,axis=0),ka,mean(Doats,axis=0),ka,mean(Doatsp,axis=0))
plot(k,Dmaxth,ka,mean(Dac,axis=0),ka,mean(Doats,axis=0))
xlabel(r'$ka$')
ylabel(r'$D$')
legend((r'$<D_\max>$ th.',r'$<D_\max>$ FAR',r'$<D_\max>$ OATS'),loc=4)
grid('on')



r=Doatsp/Dacp
rm=mean(Doatsp,axis=0)/mean(Dacp,axis=0)

r25=zeros(len(f))
r975=zeros(len(f))
for i in range(0,len(f)):
	r25[i]=prctile(r[i,:],(2.5))
	r975[i]=prctile(r[i,:],(97.5))

figure(2)
plot(f/1e6,mean(r,axis=1),'k-',f/1e6,r25,'k--',f/1e6,r975,'k--')
legend((r'$<D_\max>_{oats}/<D_\max>_{far}$',r'95 % CI'),loc=1)
xlabel(r'$f$ (MHz)')

grid('on')
show()



