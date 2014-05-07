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
f = array(arange(50e6,3e9+50e6,50e6))
np=360
nt=180
dtheta = pi/nt
dphi = (2*pi)/np
#measurement points
phi=linspace(0,2*pi,np)
theta=linspace(0,pi,nt)#arccos(2*rand(M,1)-1)



MC=500
Dac=zeros(len(f))
Doats=zeros(len(f))

#for o in range(0,MC):
#	print o
Ethac=zeros((len(theta),len(phi),len(f)),'complex')
Ephac=zeros((len(theta),len(phi),len(f)),'complex')
Erac=zeros((len(theta),len(phi),len(f)),'complex')
Ethoats=zeros((len(theta),len(phi),len(f)),'complex')
Ephoats=zeros((len(theta),len(phi),len(f)),'complex')
Eroats=zeros((len(theta),len(phi),len(f)),'complex')
TH,PH=meshgrid(theta,phi)
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
for i in range(0,len(PH[:,0])):
	for j in range(0,len(TH[0,:])):
		X=R*cos(PH[i,j])*sin(TH[i,j])
		Y=R*sin(PH[i,j])*sin(TH[i,j])
		Z=R*cos(TH[i,j])
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
		Ethac[j,i,:]= Exx*cos(TH[i,j])*cos(PH[i,j])+Eyy*cos(TH[i,j])*sin(PH[i,j])-Ezz*sin(TH[i,j])
		Ephac[j,i,:]= -Exx*sin(PH[i,j])+Eyy*cos(PH[i,j])
		Erac[j,i,:] = Exx*sin(TH[i,j])*cos(PH[i,j])+Eyy*sin(TH[i,j])*sin(PH[i,j])+Ezz*cos(TH[i,j])	

#OATS

h=1
I1=concatenate((x,y,z+h,tilt,azimut,amplitude,phas), axis=1)
I2=concatenate((x,y,-(z+h),tilt,azimut+pi,amplitude,phas), axis=1)
Ioats=vstack((I1,I2))
#Free space caracterisation (perfect anechoic chamber)
for i in range(0,len(PH[:,0])):
	for j in range(0,len(TH[0,:])):
		X=R*cos(PH[i,j])*sin(TH[i,j])
		Y=R*sin(PH[i,j])*sin(TH[i,j])
		Z=R*cos(TH[i,j])
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
		L =tile(Ioats[:,5],(len(f),1))*1/dp*ld*fp.T/c/2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
		Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
		Ezz = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
		Ethoats[j,i,:]= Exx*cos(TH[i,j])*cos(PH[i,j])+Eyy*cos(TH[i,j])*sin(PH[i,j])-Ezz*sin(TH[i,j])
		Ephoats[j,i,:]= -Exx*sin(PH[i,j])+Eyy*cos(PH[i,j])
		Eroats[j,i,:] = Exx*sin(TH[i,j])*cos(PH[i,j])+Eyy*sin(TH[i,j])*sin(PH[i,j])+Ezz*cos(TH[i,j])

for u in range(0,len(f)):
	Fa2ac=abs(Ephac[:,:,u])**2+abs(Ethac[:,:,u])**2
	Faac=Fa2ac/Fa2ac.max()
	omegaac = (Faac*sin(TH.T)*dtheta*dphi).sum()
	Dac[u] = 4*pi/omegaac
	Fa2oats=abs(Ephoats[:,:,u])**2+abs(Ethoats[:,:,u])**2
	Faoats=Fa2oats/Fa2oats.max()
	omegaoats = (Faoats*sin(TH.T)*dtheta*dphi).sum()
	Doats[u] = 4*pi/omegaoats



#savetxt('Doats.txt', Doats)
#savetxt('Dac.txt', Dac)
#savetxt('f.txt', f)



figure(1)
plot(f,Dac,f,Doats)
grid('on')


figure(2)
plot(f,Doats/Dac)
grid('on')
show()

#rrr=Doats/Dac
#figure(2)
#plot(f,mean(rrr,axis=0),f,prctile(rrr,(2.5)),f,prctile(rrr,(97.5)))
#grid('on')
#show()

#plot(f,Dac[3,:],f,Doats[3,:])
#grid('on')
#show()
#
#


#figure(1)
#for u in range(0,len(f)):
#	subplot(121)
#	im=imshow(abs(Ethac[:,:,u]))
#	cbar = colorbar(im, orientation='horizontal')
#	subplot(122)
#	im2=imshow(abs(Ephac[:,:,u]))
#	cbar2 = colorbar(im, orientation='horizontal')
#	#ticks=[-1, 0, 1]
#	show()
#	close()

#import os
#import sys#files = []##for u in range(0,len(f)):
#	close()
#	Fa2ac=abs(Ephac[:,:,u])**2+abs(Ethac[:,:,u])**2
#	Faac=Fa2ac/Fa2ac.max()
#	omega = (Fa*sin(TH.T)*dtheta*dphi).sum()
#	Dac = 4*pi/omega
#	print Dac
#	fig = figure(num=6, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
#	fig.suptitle(r'%2.2f MHz, $D=$%2.2f' %(int(round(f[u]/1e6)),Dac), fontsize=14, fontweight='bold')
#	#text(0.5, 0.5,'matplotlib',horizontalalignment='center', verticalalignment='center')
#	subplot(221)
#	im=imshow(abs(Ethac[:,:,u]))
#	xlabel(r'$\phi$ ($^\circ$)')
#	ylabel(r'$\theta$ ($^\circ$)')
#	cbar = colorbar(im, orientation='horizontal')
#	cbar.set_label(r'V/m')
#	title(r'$E_\theta$')
#	subplot(222)
#	im2=imshow(abs(Ephac[:,:,u]))
#	xlabel(r'$\phi$ ($^\circ$)')
#	ylabel(r'$\theta$ ($^\circ$)')
#	title(r'$E_\phi$')
#	cbar2 = colorbar(im, orientation='horizontal')
#	cbar2.set_label(r'V/m')
#	subplot(223)
#	hist1=hist(list(abs(Ethac[:,:,u]).flatten(1)),30)
#	xlabel(r'V/m')
#	ylabel(r'Occ.')
#	title(r'$E_\theta$')
#	#xlim(0,70)
#	subplot(224)
#	hist2=hist(list(abs(Ephac[:,:,u]).flatten(1)),30)
#	title(r'$E_\phi$')
#	#xlim(0,70)
#	xlabel(r'V/m')
#	ylabel(r'Occ.')
#	fname = '_tmp%03d.png'%u#	print 'Saving frame', fname#	fig.savefig(fname)#	files.append(fname)
#	close()
##