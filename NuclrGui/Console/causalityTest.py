from pylab import *
from scipy.misc.common import derivative

mpi = 135.0
mn = 939.0

Eb = -15.8
Asymm = 32.0
dim = mpi**3
n0 = 0.5

def U1(n):
    return n

def U2(n):
    return n/(1 - 0.02*n*exp(-(6.0/n)**4.0))

def E(n, u):
    return dim*n*(mn + Eb*u(n/n0) * (2.2 - u(n/n0))/(1 + 0.2*u(n/n0)))

def P(n, u):
    return n*derivative(lambda z: E(z, u), n, order = 5, dx=1e-5) - E(n, u)

def C(n, u):
    num = derivative(lambda z: P(z, u), n, order = 5)
    den = derivative(lambda z: E(z, u), n, order = 5)
    return num/den

print C(8*n0, U2)