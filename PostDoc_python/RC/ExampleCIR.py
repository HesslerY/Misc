#CIR.py
from numpy import *
from pylab import *
from scipy import fftpack

#calculation of the images
Tm = 1e-6
l = 8.7
p = 3.7
h = 2.9

dims = array([l,p,h])
c = 299792458

dmax = Tm*c
Tp = Tm+3./c*dims.max()

X = 1.
Y = 2.
Z = 1.

tilt = pi/2-math.acos(sqrt(2./3));
azimut = pi/4;

order = round(dmax/dims.min())

Memo = pi/6*(Tp*c)**3/l/p/h*6*8

POSp = array([X,Y,Z, 0, tilt, azimut]);
#1D

for i in range(1,int(order)):
	POSp = vstack([POSp,[2*i*l-X, Y, Z, abs(2*i)-1, pi-tilt, 2*pi-azimut],[2*i*l+X, Y, Z, abs(2*i), tilt, azimut]])
	#POSp=append(POSp,[[2*i*l-X, Y, Z, abs(2*i)-1, pi-tilt, 2*pi-azimut],[2*i*l+X, Y, Z, abs(2*i), tilt, azimut]],axis=0)

POSi=POSp.copy();
POSi[:,1] = array([p*ones((1,len(POSp)))-POSi[:,1]])
POSi[:,4] = array([pi*ones((1,len(POSp)))-POSi[:,4]])
POSi[:,5] = array([pi*ones((1,len(POSp)))-POSi[:,5]])

#2D
#POSP= empty((1,6))
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
	#U = where(POS2[:,6]<c*Tp)
	#if len(U)>2:
	#POS2 = POS2[0:max(max(U)),:]
	V = where(POS2[:,6]>c*Tp)
	POS2=delete(POS2, s_[min(min(V)):], axis=0)
	POS2 = delete(POS2, s_[6], axis=1)
	POSP = concatenate((POSP,POS2),axis=0)
	

POSI = POSP.copy();
POSI[:,2] = array([h*ones((1,len(POSI)))-POSI[:,2]])
POSI[:,5] = mod(array([pi*ones((1,len(POSI)))+POSI[:,5]]),2*pi)

#3D

POSPP = POSP.copy()
#POSPP[:,2] = array([k*h*ones((1,len(POSP)))]+POSP[:,2])
#POSPP[:,3] = array([abs(k)*ones((1,len(POSP)))]+POSP[:,3])
POSII=POSI.copy()
POSII[:,2] = array([h*ones((1,len(POSI)))]+POSI[:,2])
POSII[:,3] = array([ones((1,len(POSI)))]+POSI[:,3])
POS3 = POSPP.copy()
POS3 = concatenate((POS3,POSII),axis=0)
dist = (POS3[:,0]**2+POS3[:,1]**2+POS3[:,2]**2)**.5
POS3 = concatenate((POS3,dist.reshape(-1,1)),axis=1)
POS3 = POS3[POS3[:,6].argsort(),]
W = where(POS3[:,6]>c*Tp)
POS3=delete(POS3, s_[min(min(W)):], axis=0)#POS = POS3[0:max(max(U)),:]
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

#Example.py

#Channel impulse response calculation
R=0.998

#measurement point
X_1=array([4.5])
Y_1=array([3])
Z_1=array([1])

def CIR(a,b,c):
	#1/8
	Sx1,Sy1,Sz1 = CIR8th(a,b,c)
	2/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx2,Sy2,Sz2 = CIR8th(a,b,c)
	POS[:,2]=-POS[:,2]
 	#3/8
	POS[:,1]=-POS[:,1]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	Sx3,Sy3,Sz3 = CIR8th(a,b,c)
	#4/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx4,Sy4,Sz4 = CIR8th(a,b,c)
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	#5/8
	POS[:,3]=array([POS[:,3]-ones((1,len(POS)))])
	POS[:,0]=-POS[:,0]
	Sx5,Sy5,Sz5 = CIR8th(a,b,c)
	#6/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx6,Sy6,Sz6 = CIR8th(a,b,c)
	POS[:,2]=-POS[:,2]
	POS[:,4]=array([pi*ones((1,len(POS)))-POS[:,4]])
	#7/8
	POS[:,1]=-POS[:,1]
	POS[:,5]=-POS[:,5]
	Sx7,Sy7,Sz7 = CIR8th(a,b,c)
	#8/8
	POS[:,3]=array([ones((1,len(POS)))+POS[:,3]])
	POS[:,2]=-POS[:,2]
	POS[:,5]=array([pi*ones((1,len(POS)))+POS[:,5]])
	Sx8,Sy8,Sz8 = CIR8th(a,b,c)
	POS[:,0]=-POS[:,0]
	POS[:,1]=-POS[:,1]
	POS[:,2]=-POS[:,2]
	POS[:,3]=array([-3.*ones((1,len(POS)))+POS[:,3]])
	POS[:,5]=mod(POS[:,5],2*pi)
	Sx=Sx1+Sx2+Sx3+Sx4+Sx5+Sx6+Sx7+Sx8
	Sy=Sy1+Sy2+Sy3+Sy4+Sy5+Sy6+Sy7+Sy8
	Sz=Sz1+Sz2+Sz3+Sz4+Sz5+Sz6+Sz7+Sz8
	return Sx,Sy,Sz




def CIR8th(q,s,w):
	Sx8th=zeros((len(q),len(ts)))
	Sy8th=zeros((len(q),len(ts)))
	Sz8th=zeros((len(q),len(ts)))
	for i in range(0,len(q)):
		for j in range(0,len(POS)):
			DX = q[i]-POS[j,0]
			DY = s[i]-POS[j,1]
			DZ = w[i]-POS[j,2]
			r = sqrt(DX**2+DY**2+DZ**2)
			u=int(r/c/Tm*N)
			if u<N:
				ca    = cos(POS[j,4])
				sa    = sin(POS[j,4])
				cb    = cos(POS[j,5])
				sb    = sin(POS[j,5])
				rx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ
				ry = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
				rz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
				DXY=sqrt(DX**2+DY**2)
				rxy = sqrt(rx**2+ry**2)
				costheta = rz/r
				sintheta = rxy/r
				cosphi   = rx/rxy
				sinphi   = ry/rxy
				L =R**(POS[j,3])*1/r      
				Sx8th[i,u] = Sx8th[i,u]+L*(((-sb)**2+(1-(-sb)**2)*ca)*(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))*(-sintheta*costheta*sinphi)+(-cb*sa)*(-sintheta*(-sintheta)));
				Sy8th[i,u] = Sy8th[i,u]+L*((-sb*cb*(1-ca))*(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)*(-sintheta*costheta*sinphi)+(-sb*sa)*(-sintheta*(-sintheta)));
				Sz8th[i,u] = Sz8th[i,u]+L*((cb*sa)*(-sintheta*costheta*cosphi)+(sb*sa)*(-sintheta*costheta*sinphi)+ca*(-sintheta*(-sintheta)));
	return Sx8th,Sy8th,Sz8th

#sampling
f0=1e9
N=int(3*Tm*f0)
ts=Tm/N*arange(1,N+1)
Sx,Sy,Sz=CIR(X_1,Y_1,Z_1)


for i in range(0,len(Sz[:,0])):
	figure(i)
	subplot(311)
	l = plot(ts/1e-6,(Sx[i,:]))
	grid(True)
	title('Channel impulse response')
	ylabel('$E_x$')	
	subplot(312)
	l = plot(ts/1e-6,(Sy[i,:]))
	grid(True)
	ylabel('$E_y$')
	subplot(313)
	l = plot(ts/1e-6,(Sz[i,:]))
	ylabel('$E_z$')
	grid(True)
	xlabel('time ($\mu$s)')

def nextpow2(v):
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v + 1

Fs = N/Tm #sampling frequency
periode = 1/Fs                     # sample length
L = N                  #Number of points
tt = [0:L-1]*periode                # time...
NFFT = 2^nextpow2(N) # Next power of 2 from length of y

FFTx=zeros((len(Sz[:,0]),NFFT/2-1),'complex')
FFTy=zeros((len(Sz[:,0]),NFFT/2-1),'complex')
FFTz=zeros((len(Sz[:,0]),NFFT/2-1),'complex')

for i in range(0,len(Sz[:,0])):
	Yx = fft(Sx[i,:],NFFT)
	Yy = fft(Sy[i,:],NFFT)
	Yz = fft(Sz[i,:],NFFT)
	f = Fs/2*linspace(0,1,NFFT/2)
	FFTx[i,:] = Yx[1:NFFT/2]
	FFTy[i,:] = Yy[1:NFFT/2]
	FFTz[i,:] = Yz[1:NFFT/2]
	freq = f

for i in range(0,len(Sz[:,0])):
	figure(i+2)
	subplot(311)
	l = plot(freq[0:len(FFTx[i,:])]/1e6,20*log(abs(FFTx[i,:])))
	grid(True)
	title('Frequency response')
	ylabel('$FFT_x$')	
	subplot(312)
	l = plot(freq[0:len(FFTy[i,:])]/1e6,20*log(abs(FFTy[i,:])))
	grid(True)
	ylabel('$FFT_y$')
	subplot(313)
	l = plot(freq[0:len(FFTz[i,:])]/1e6,20*log(abs(FFTz[i,:])))
	ylabel('$FFT_z$')
	grid(True)
	xlabel('frequency (MHz)')

show()
