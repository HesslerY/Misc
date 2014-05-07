#! /Library/Frameworks/Python.framework/Versions/2.7/bin

#Brownian.py
from numpy import *
from numpy.random import *
from pylab import *

N=10000
M=100000

r=rand(N,1)
phi=2*pi*rand(N,1)
cp=r*cos(phi)
sp=r*sin(phi)

X=cp.cumsum(axis=0)
Y=sp.cumsum(axis=0)

plot(X,Y)
grid(True)
xlabel('$x$')
ylabel('$y$')
show()
