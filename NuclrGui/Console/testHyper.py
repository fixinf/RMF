import eosCore as eos
from RMF.Wrapper import Wrapper
from pylab import *

HyperC = eos.Walecka(sqrt(164.462), sqrt(54.6), sqrt(121.69), 0.0202, 0.047, 0)

print eos.test_l0(2.0, HyperC)
print eos.test_sm(2.0, HyperC)

x = arange(0.1, 3.0, 0.01)
def n_l(n, C):
    np = eos.np_eq(n, C)
    f = eos.f_eq(n - np, np, C)
    mu_e = eos.mu_e(n, np, f, C)
    res = 0
    if mu_e > 0.5:
        res += (mu_e**2 - 0.5**2)**(1.5) / (3*pi**2)
    if mu_e > 105.0:
        res += (mu_e**2 - 105**2)**(1.5) / (3*pi**2)
    return res/134.92**3
        
print eos.test_xm(5.0, HyperC)