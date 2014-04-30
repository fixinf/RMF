import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max
print 'hey!'
C0 = eos.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
C = eos.gc_OR5_S(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

wr0 = Wrapper(C0, './tab/KVOR.dat')
wr0.Solve()
wr0.tabulate()

C.set_f0(0.21)
C.set_gp_beta(2.8)
C.set_gp_d(0.0)
C.set_gp_e(0.0)
C.set_gp_gamma(2.0)
C.set_gp_z(0.65)

wr = Wrapper(C, './tab/OR5_S.dat')
wr.Solve()
wr.tabulate()
M0, R0 = wr0.star(0.5, 5.0, 0.25)
M, R = wr.star(0.5, 5.0, 0.25)
print max(M0), max(M)
plot(R0, M0, R, M)
show()


