import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max

C = eos.gc_OR5_S(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

C.set_f0(0.21)
C.set_gp_beta(2.8)
C.set_gp_d(0.0)
C.set_gp_e(0.0)
C.set_gp_z(0.8)
wr = Wrapper(C, './tab/rho_test.dat')
glist = [2.0, 3.0, 4.0]
n = arange(0.01, 8*wr.n0, wr.n0/20)
fig, ax = subplots(1, 4) 
mmax = []
for gamma in glist:
    C.set_gp_gamma(gamma)
    wr.filename = './tab/rho_test/gamma=%.1f' % gamma
    print wr.filename
    wr.Solve()
    wr.tabulate()
    np_eq = map(lambda z: wr.np_eq(z)/z, n)
    ax[0].plot(n, np_eq)
    f_eq = map(lambda z: wr.f_eq(z - wr.np_eq(z), wr.np_eq(z)), n)
    ax[1].plot(n, f_eq)
    p = wr.P(n)
    ax[2].plot(n, p)
    M, R = wr.star(wr.n0, wr.n0*10, wr.n0/5)

    mmax.append(max(M))
print mmax, gamma
ax[3].plot(array(glist), array(mmax))    
show()

# C.set_gp_gamma(2.0)
# blist = [2.8, 3.5, 4.0, 4.5]
# fig, ax = subplots(1, 2) 
# 
# for beta in blist:
#     C.set_gp_beta(beta)
#     wr.Solve()
#     np_eq = map(lambda z: wr.np_eq(z)/z, n)
#     ax[0].plot(n, np_eq)
#     f_eq = map(lambda z: wr.f_eq(z - wr.np_eq(z), wr.np_eq(z)), n)
#     ax[1].plot(n, f_eq)
# 
# show()
# raw_input('press a key')