from causalityTest import *
import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from scipy import optimize

Cs =  eos.gc_OR5_S_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj2 = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj3 = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

KV =  eos.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
KV.set_f0(0.195)
wrKV = Wrapper(KV, '')
wrKV.Solve()

f0 = 0.2

Csj.set_gp_alpha(0.0)
Csj.set_gp_z(0.65)
Csj.set_gp_gamma(2.0)
Csj.set_gp_beta(2.8)
Csj.set_gp_w(0.0)
Csj.set_gp_e(-0.0)
Csj.set_f_tilde(10.0)
#--------------------------------
acj = 0.0
Csj.set_gp_d(0.0)

Csj.set_omega_a(0.0)
Csj.set_fcut_omega(0.6)
Csj.set_omega_c(0.0)

Csj.set_gp_a(0)
Csj.set_gp_fcut(0.0)

#---------------------------------

Csj.set_omega_b(0.0)

Csj.set_rho_a(0.0)
Csj.set_rho_b(0.0)

#------------------------------
Csj.set_f0(0.2)
wrj = Wrapper(Csj,'')
wrj.Solve(-15.8, 250.0, 28.0)

print wrj.E_symm(wrj.n0)
print E(n0, U1)/(dim*n0)

n = arange(0.001, 5.0, 0.01)

# plot(n/wrj.n0, wrj.E_symm(n))
# plot(n/n0, map(lambda z: E(z, U1)/(dim*z) - mn, n))
# show()

const = 197.33*(0.16/n0)**(4.0/3.0)*mpi**(-4.0)

UpperY = [15.639, 60, 140, 210]
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]
n_Constr = arange(0.5, 4.5, 0.01)
# semilogy(UpperX, UpperY, LowerX, LowerY, color='red', linewidth=2.0)
# semilogy(n_Constr/wrj.n0, wrj.P_symm(n_Constr))
# semilogy(n_Constr/n0, map(lambda z:const*P(z, U2), n_Constr))
# semilogy(n_Constr/n0, map(lambda z:const*P(z, U1), n_Constr))
# 
# ylim(1.0, 400.0)
# xlim(1, 6.25)
# show()

n_fit = linspace(6.0, 8.0, 5)
print n_fit/wrj.n0
def set_params(wr, f_o=0.4, a_o=0.0, f_phi=0.4, a_phi=0.0, a_c=0.0, d=0.0):
#     wr.C.set_gp_alpha(alpha)
    wr.C.set_fcut_omega(f_o)
    wr.C.set_omega_a(a_o)
    wr.C.set_gp_d(d)
    wr.C.set_gp_a(a_phi)
    wr.C.set_gp_fcut(f_phi)
    wr.Solve(-15.8, 250.0, 28.0)

def lsq_fun(x):
    print x
    set_params(wrj, a_o = x[0], a_phi=x[1], a_c=x[2], d = x[3])
#     semilogy(n_Constr/wrj.n0, wrj.P_symm(n_Constr))
#     semilogy(n_Constr/n0, map(lambda z:const*P(z, U2), n_Constr))
#     show()
    print array(wrj.P_symm(n_fit*n0))
    print array(map(lambda z: const*P(z, U2),n0*n_fit))
    res = array(wrj.P_symm(n_fit*n0)) - array(map(lambda z: const*P(z, U2), n_fit*n0))
    print res
    return res


x = optimize.leastsq(lsq_fun, [10.0, -1.0, -300, -1.0])[0]
print x
set_params(wrj, a_o=x[0], a_phi=x[1])
# set_params(wrj, a_o=0, a_phi=-100)
semilogy(UpperX, UpperY, LowerX, LowerY, color='red', linewidth=2.0)
semilogy(n_Constr/wrj.n0, wrj.P_symm(n_Constr))
set_params(wrj)
semilogy(n_Constr/wrj.n0, wrj.P_symm(n_Constr))
semilogy(n_Constr/n0, map(lambda z:const*P(z, U2), n_Constr))
xlim(1, 9)
ylim(1.0, 1200)
show()


# print wrj.C.C_s**2, wrj.C.C_o**2, wrj.C.C_r**2,wrj.C.b,wrj.C.c



