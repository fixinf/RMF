import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max
from gtk._gtk import Orientation
from scipy import optimize
Csj = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

acj = 0.5
Csj.set_gp_alpha(1.0)
Csj.set_gp_z(0.65)
Csj.set_gp_gamma(2.0)
Csj.set_gp_beta(2.8)
Csj.set_gp_w(1.0)
Csj.set_gp_e(-0.0)
Csj.set_gp_d(1.0)
Csj.set_f_tilde(10.0)

Csj.set_gp_a(0.0)
Csj.set_gp_fcut(0.61)

Csj.set_omega_a(0)
Csj.set_omega_b(0.0)
Csj.set_fcut_omega(0.6)
Csj.set_rho_a(0.0)
Csj.set_rho_b(10.0)

wrj = Wrapper(Csj, 'Csj.tab')
n0 = wrj.n0

Csj.set_f0(0.24)
Csj.set_gp_fcut(0.35)
Csj.set_fcut_omega(0.55)

def set_params(x, wr):
    alpha = x[0]
    d = x[1]
    a_phi = x[2]
    a_omega = x[3]
#     f0 = x[4]
    wr.C.set_omega_a(a_omega)
    wr.C.set_gp_d(d)
    wr.C.set_gp_a(a_phi)
#     wr.C.set_f0(f0)
    for a in linspace(wr.C.gp_alpha, alpha, 3):
        wr.C.set_gp_alpha(a)
        wr.Solve()
        

UpperY = array([15.639, 60, 140, 210])
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= array([2, 2.5, 3.5, 4, 4.5])
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]


# wrj.Solve()
# print wrj.star_notab(0.5, 4, 0.5)

def lsq_fun(x):
    print x
    set_params(x, wrj)
    res = array(UpperY) - array(wrj.P_symm(array(UpperX)*n0))
#     wrj.tabulate()
#     M,R = wrj.star(2*n0, 8*n0, n0/10)
#     mmax = max(M)
    return res
#     return append(res, array([mmax - 2.2]))
# 
# print lsq_fun([1.0, 0.0, 0.0])
# print lsq_fun([1.1, 0.1])
print '--------------------------------------------------------------'
x = optimize.leastsq(lsq_fun, [1.4, -5.5, 2.5, -4.5], full_output=True, factor=0.2, xtol=1e-3)[0]
#x = optimize.anneal(lsq_fun, [1.3, 1.0, 0.0, 0.0, 0.24], full_output=True, lower=[0, -3.0, -2.0, -2.0, 0.21], 
#                     upper=[2.0, 2.0, 2.0, 2.0, 0.3])[0]
# x = optimize.fmin_cg(lsq_fun, asarray((1.3, 1.0, 0.0, 0.0)), full_output=True)[0]

print x

nConstr = arange(n0, 5*n0, n0/10)
set_params(x, wrj)
semilogy(UpperX, UpperY, LowerX, LowerY,color='red')
semilogy(nConstr/n0, wrj.P_symm(nConstr))
ylim([1.0, 400.0])
show()

wrj.tabulate()
nstar = arange(2*n0, 8*n0, n0/10)
M,R = wrj.star(2*n0, 8*n0, n0/10)
print x
print max(M)



