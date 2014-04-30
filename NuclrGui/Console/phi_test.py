import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max

C = eos.gc_OR5_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
C0 = eos.gc_OR5_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Cs =  eos.gc_OR5_S_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

KV =  eos.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
KV.set_f0(0.21)
wrKV = Wrapper(KV, '')
wrKV.Solve()

f0 = 0.24

C.set_f0(f0)
C.set_gp_beta(2.8)
C.set_gp_z(0.65)
C.set_gp_fcut(0.45)
C.set_gp_a(1.0)
C.set_gp_alpha(1.0)
C.set_gp_gamma(2.0)
C.set_gp_w(1.0)

C0.set_f0(f0)
C0.set_gp_beta(2.8)
C0.set_gp_z(0.65)
C0.set_gp_fcut(0.4)
C0.set_gp_a(0.0)
C0.set_gp_alpha(1.0)
C0.set_gp_w(1.0)
C0.set_gp_gamma(2.0)

Cs.set_f0(f0)
Cs.set_gp_beta(2.8)
Cs.set_gp_z(0.65)

Cs.set_gp_alpha(1.0)
Cs.set_gp_gamma(2.0)
Cs.set_gp_w(1.0)

Cs.set_gp_e(-8)
Cs.set_gp_d(-15)
Cs.set_gp_a(1.0)
Cs.set_gp_fcut(0.22)
acs = 1.8


Csj.set_f0(0.27)

acj = 1.2
Csj.set_gp_alpha(1.0)
Csj.set_gp_z(0.65)


Csj.set_gp_gamma(2.0)
Csj.set_gp_beta(2.8)

Csj.set_gp_w(1.0)

Csj.set_gp_e(0.0)
Csj.set_gp_d(-4.5)
Csj.set_f_tilde(0.6)

Csj.set_gp_a(1.5)
Csj.set_gp_fcut(0.35)

Csj.set_omega_a(-3.0)
Csj.set_omega_b(0.0)
Csj.set_fcut_omega(0.55)

Csj.set_rho_a(0.0)
Csj.set_rho_b(10.0)

f_tilde_list = [1.5]

wr = Wrapper(C, './tab/rho_test.dat')
wr0 = Wrapper(C0, './tab/rho_test.dat')

wrs = Wrapper(Cs, '')
wrj = Wrapper(Csj,'')

wrlist = [wrs, wrj, wrKV]

x = arange(0.0, 1.0, 0.01)
n = arange(0.01, 8*wr.n0, wr.n0/10)

fig, ax = subplots(3,3)

a_seq = arange(1.0, 1.2, 0.05)
for a in a_seq:
    for f in linspace(0.195, f0, 5):
        print 'a=',a,'f_0=',f0
        for wr in wrlist:
            if wr != wrKV:
                wr.C.set_gp_alpha(a)
                wr.C.set_f0(f)
                wr.Solve()

for a in linspace(1.1, acs, 3):
    print 'acs=',a
    wrs.C.set_gp_alpha(a)
    wrs.Solve()
  
for a in linspace(1.1, acj,3):
    print 'acj=', a
    wrj.C.set_gp_alpha(a)
    wrj.Solve()
    
k = wr.n0/0.16    
UpperY = [15.639, 60, 150, 200]
UpperX = [k*0.32, k*0.43, k*0.574, k*0.68]
LowerX= [k*0.32, k*0.4, k*0.64, k*0.68 ]
LowerY = [7.0, 17.49, 52.0, 52.0]
ax[1][1].semilogy(UpperX, UpperY, LowerX, LowerY)

for wr in wrlist:
    if wr == wrj:
        for f_tilde in f_tilde_list:
            wr.C.set_f_tilde(f_tilde)
            wr.Solve()
            ax[0][0].plot(x, map(lambda z: wr.C.phi_n(z),x))
            ax[0][1].plot(n, wr.P(n))
            ax[0][2].plot(n, map(lambda z: wr.f_eq(z - wr.np_eq(z),wr.np_eq(z)), n))
            ax[1][0].plot(x, map(lambda z: wr.C.eta_o(z),x))  
            ax[1][1].semilogy(n, wr.P_symm(n))
            ax[1][2].plot(x, map(lambda z: wr.C.eta_r(z),x))
            ax[2][0].plot(n, wr.E_symm(n)) 
  
    else:
            ax[0][0].plot(x, map(lambda z: wr.C.phi_n(z),x))
            ax[0][1].plot(n, wr.P(n))
            ax[0][2].plot(n, map(lambda z: wr.f_eq(z - wr.np_eq(z),wr.np_eq(z)), n))
            ax[1][0].plot(x, map(lambda z: wr.C.eta_o(z),x))
            ax[1][1].semilogy(n, wr.P_symm(n))
            ax[1][2].plot(x, map(lambda z: wr.C.eta_r(z),x)) 
            ax[2][0].plot(n, wr.E_symm(n)) 

show()


mmax = []
              
for i, wr in enumerate(wrlist):
    wr.filename='%i.dat'%i
    wr.tabulate()
    M, R = wr.star(wr.n0, 10*wr.n0, wr.n0/10)
    mmax.append(max(M))
print mmax
#     C.set_gp_gamma(gamma)
#     wr.filename = './tab/rho_test/gamma=%.1f' % gamma
#     print wr.filename
#     wr.Solve()
#     wr.tabulate()
#     np_eq = map(lambda z: wr.np_eq(z)/z, n)
#     ax[0].plot(n, np_eq)
#     f_eq = map(lambda z: wr.f_eq(z - wr.np_eq(z), wr.np_eq(z)), n)
#     ax[1].plot(n, f_eq)
#     p = wr.P(n)
#     ax[2].plot(n, p)
#     M, R = wr.star(wr.n0, wr.n0*10, wr.n0/5)
# 
#     mmax.append(max(M))
# print mmax, gamma
# ax[3].plot(array(glist), array(mmax))    
# show()