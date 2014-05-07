#PT.py
from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams

f = genfromtxt('f_500_n30_R25cm.txt', unpack=True)
Dac=genfromtxt('Dac_500_n30_R25cm.txt', unpack=True)
Doats=genfromtxt('Doats_500_n30_R25cm.txt', unpack=True)

c = 299792458
ka=2*pi*f/c*.25
k=array(arange(1,70,1))
Dmaxth=.5*(0.577+log(4*k**2+8*k)+1/(8*k**2+16*k))

figure(1)
plot(ka,mean(Dac,axis=1),'k--',ka,mean(Doats,axis=1),'k*-',k,Dmaxth,'r-')
xlim((0,70))
xlabel(r'$ka$')
ylabel(r'$D$')
legend((r'$<D_\max>$ FAR',r'$<D_\max>$ OATS',r'$<D_\max>$ th.'),loc=4)
grid('on')
show()

r=Doats/Dac

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