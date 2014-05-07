#!/usr/bin/env python

from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams
import time
import os
import sys
files = []

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
c = 299792458
R = 10.
f = array(arange(10e6,3e9+10e6,10e6))
#f = array([1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])*1e6
np=360
nt=180
dtheta = pi/nt
dphi = (2*pi)/np
#measurement points
phi=linspace(0,2*pi,np)
theta=linspace(0,pi,nt)#arccos(2*rand(M,1)-1)
Dac=zeros(len(f))
#MC=500
#Dac=zeros((MC,len(f)))
#Doats=zeros((MC,len(f)))
start = time.time()
#for o in range(0,MC):  
Ethac=zeros((len(theta),len(phi),len(f)),'complex')
Ephac=zeros((len(theta),len(phi),len(f)),'complex')
#Erac=zeros((len(theta),len(phi),len(f)),'complex')
#Real OATS
H2=4.
H1=1.
h=1.
z=linspace(H1,H2,round(4*nt/pi*(arctan(H2/R)-arctan(H1/R))))
Ezoats=zeros((len(z),len(phi),len(f)),'complex')
Ethoats=zeros((len(z),len(phi),len(f)),'complex')
#Eroats=zeros((len(z),len(phi),len(f)),'complex')
#FAR
TH,PH=meshgrid(theta,phi)
R_eut=.5 #m
n=20#number of dipoles
I=zeros((n,7))
theta_eut=arccos(2*rand(n,1)-1)
phi_eut=2*pi*rand(n,1)
xx=R_eut*cos(phi_eut)*sin(theta_eut)
yy=R_eut*sin(phi_eut)*sin(theta_eut)
zz=R_eut*cos(theta_eut)
tilt=arccos(2*rand(n,1)-1)
azimut=2*pi*rand(n,1)
ld=.1      
amplitude=rand(n,1)
phas=2*pi*rand(n,1)
I=concatenate((xx,yy,zz,tilt,azimut,amplitude,phas), axis=1)
I1=concatenate((xx,yy,zz+h,tilt,azimut,amplitude,phas), axis=1)
I2=concatenate((xx,yy,-(zz+h),tilt,azimut+pi,amplitude,phas), axis=1)
Ioats=vstack((I1,I2))

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
                #Erac[j,i,:] = Exx*sin(TH[i,j])*cos(PH[i,j])+Eyy*sin(TH[i,j])*sin(PH[i,j])+Ezz*cos(TH[i,j])     
        for l in range(0,len(z)):
                X=R*cos(phi[i])
                Y=R*sin(phi[i])
                Z=z[l]
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
                k=len(z)-l-1
                L =tile(Ioats[:,5],(len(f),1))*1/dp*ld*(fp.T/c)**2*377#377/4/pi*2*pi*f/c*repmat(I(:,6),1,length(f))*ld/dp; %Amplitude & free space attenuation
                Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
                Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
                Ezoats[k,i,:] = sum(exp(1j*phase)*L*tile((((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
                Ethoats[k,i,:] = -sin(phi[i])*Exx+cos(phi[i])*Eyy
                #Eroats[k,i,:] = cos(phi[i])*Exx+sin(phi[i])*Eyy

for u in range(0,len(f)):
        close()
        pac=abs(Ephac[:,:,u])**2+abs(Ethac[:,:,u])**2
        rPac=pac.max()/pac.mean()
        poats=abs(Ezoats[:,:,u])**2+abs(Ethoats[:,:,u])**2
        rPoats=poats.max()/poats.mean()
        rE=vstack((abs(Ezoats[:,:,u]),abs(Ethoats[:,:,u]))).mean()/vstack((abs(Ephac[:,:,u]),abs(Ethac[:,:,u]))).mean()
        fig = figure(num=6, figsize=(10, 7), dpi=200, facecolor='w', edgecolor='k')
        fig.suptitle(r'%2.2f MHz, $\max(P_{r\textrm{ac}})/<P_{r\textrm{ac}}>=$%2.2f, $\max(P_{r\textrm{oats}})/<P_{r\textrm{oats}}>=$%2.2f, $<E_{r\textrm{oats}}>/<E_{r\textrm{ac}}>=$%2.2f' %(int(round(f[u]/1e6)),rPac,rPoats,rE), fontsize=14, fontweight='bold')
        #text(0.5, 0.5,'matplotlib',horizontalalignment='center', verticalalignment='center')
        subplot(221)
        im=imshow(abs(Ethac[:,:,u]))
        xlabel(r'$\phi$ ($^\circ$)')
        ylabel(r'$\theta$ ($^\circ$)')
        cbar = colorbar(im, orientation='horizontal')
        cbar.set_label(r'V/m')
        title(r'$E_\theta$')
        subplot(222)
        im2=imshow(abs(Ephac[:,:,u]))
        xlabel(r'$\phi$ ($^\circ$)')
        ylabel(r'$\theta$ ($^\circ$)')
        title(r'$E_\phi$')
        cbar2 = colorbar(im2, orientation='horizontal')
        cbar2.set_label(r'V/m')
        subplot(223)
        im3=imshow(abs(Ezoats[:,:,u]),aspect=2)
        yticks(arange(0,len(z),20), ('4', '3', '2', '1') )
        ylabel(r'$z$ (m)')
        xlabel(r'$\theta$ ($^\circ$)')
        cbar = colorbar(im3, orientation='horizontal')
        cbar.set_label(r'V/m')
        title(r'$E_z$')
        subplot(224)
        im4=imshow(abs(Ethoats[:,:,u]),aspect=2)
        ylabel(r'$z$ (m)')
        yticks(arange(0,len(z),20), ('4', '3', '2', '1') )
        xlabel(r'$\theta$ ($^\circ$)')
        title(r'$E_\theta$')
        cbar2 = colorbar(im4, orientation='horizontal')
        cbar2.set_label(r'V/m')
        fname = '20elem%03d.png'%u
        print 'Saving frame', fname
        fig.savefig(fname)
        files.append(fname)
        close()




