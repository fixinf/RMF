import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max
from gtk._gtk import Orientation

C = eos.gc_OR5_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
C0 = eos.gc_OR5_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Cs =  eos.gc_OR5_S_F(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj2 = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
Csj3 = eos.gc_OR5_S_F_J(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)

KV =  eos.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
KV.set_f0(0.195)
wrKV = Wrapper(KV, '')
wrKV.Solve()

f0 = 0.2
Cs.set_f0(f0)
Cs.set_gp_beta(2.8)
Cs.set_gp_z(0.65)

Cs.set_gp_alpha(0.0)
Cs.set_gp_gamma(2.0)
Cs.set_gp_w(0.0)

Cs.set_gp_e(0)
Cs.set_gp_d(0)
Cs.set_gp_a(0.0)
Cs.set_gp_fcut(0.22)
acs = 0.0

Csj.set_gp_alpha(0.0)
Csj.set_gp_z(0.65)
Csj.set_gp_gamma(2.0)
Csj.set_gp_beta(2.8)
Csj.set_gp_w(1.0)
Csj.set_gp_e(0.0)

Csj.set_f_tilde(10.0)

#--------------------------------
acj = 0.0
Csj.set_gp_d(0.0)

Csj.set_omega_a(30.0)
Csj.set_fcut_omega(0.4)
Csj.set_omega_c(-330.0)


Csj.set_gp_a(10.0)
Csj.set_gp_fcut(0.4)


#---------------------------------

Csj.set_omega_b(0.0)

Csj.set_rho_a(0.0)
Csj.set_rho_b(0.0)

#------------------------------

masses = 0

f_tilde_list = [1.5]
alist = [acj]
wr = Wrapper(C, './tab/rho_test.dat')
wr0 = Wrapper(C0, './tab/rho_test.dat')

wrs = Wrapper(Cs, '')
wrj = Wrapper(Csj,'')
wrj2 = Wrapper(Csj2,'')
wrs.name='KVR'
wrKV.name = 'KVR'
wrj.name = 'KVR Modified'
wrj2.name = 'II'
wrlist = [wrj, wrs]

x = arange(0.0, 1.0, 0.01)
n = arange(0.01, 8*wr.n0, wr.n0/10)

ion()
fig, ax = subplots(3,3)

fig.set_dpi(100)

a_seq = arange(1.0, 1.2, 0.05)

for f in linspace(0.195, f0, 5):
    for wr in wrlist:
        if wr != wrKV:
            wr.C.set_f0(f)
            wr.Solve()

wrs.C.set_gp_a(0.0)
wrs.Solve(-15.8, 250.0, 28.0)  
print wrs.E_symm(wrs.n0)

for a in linspace(1.0, acj,3):
    print 'acj=', a
    wrj.C.set_gp_alpha(a)
    wrj.Solve(-15.8, 250.0, 28.0)

k = 1.0/0.16    
UpperY = [15.639, 60, 140, 210]
UpperX = [2, 2.75, 3.5, 4.5]
LowerX= [2, 2.5, 3.5, 4, 4.5]
LowerY = [7.5, 17.49, 40.0, 60.0, 60.0]
n_Constr = arange(wr.n0, 5*wr.n0, 0.01)
ax[1][1].semilogy(UpperX, UpperY, LowerX, LowerY, color='red', linewidth=2.5)
for axe in ax.flatten():
    axe.tick_params(axis='both', which='major', labelsize=8)
    
lines = []
for wr in wrlist:
    
    line, = ax[0][0].plot(x, map(lambda z: wr.C.phi_n(z),x))
    lines.append(line)
    ax[0][0].set_xlabel(r'$f$', fontsize = 18)
    ax[0][0].set_ylabel(r'$\Phi(f)$', fontsize=18, rotation='horizontal' )
    
    ax[0][1].plot(n/wr.n0, wr.P(n))
    ax[0][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[0][1].set_ylabel(r'$P_{NSM}(n)$',fontsize=18, rotation='horizontal')
    
    ax[0][2].plot(n/wr.n0, map(lambda z: wr.C.phi_n(wr.f_eq(z - wr.np_eq(z),wr.np_eq(z))), n))
    ax[0][2].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[0][2].set_ylabel(r'$\frac{m^\ast}{m_n}$',fontsize=18, rotation='horizontal')
    
    ax[1][0].plot(x, map(lambda z: wr.C.eta_o(z),x))
    ax[1][0].set_xlabel(r'$f$',fontsize=18, rotation='horizontal')
    ax[1][0].set_ylabel(r'$\eta_\omega(f)$',fontsize=18, rotation='horizontal') 
    
    ax[1][1].semilogy(n_Constr/wr.n0, wr.P_symm(n_Constr))
    ax[1][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[1][1].set_ylabel(r'$P_{SNM}(n)$',fontsize=18, rotation='horizontal')
    ax[1][1].set_ylim([1.0, 400.0])
    
    ax[1][2].plot(x, map(lambda z: wr.C.eta_r(z),x))
    ax[1][2].set_xlabel(r'$f$',fontsize=18, rotation='horizontal')
    ax[1][2].set_ylabel(r'$\eta_\rho(f)$',fontsize=18, rotation='horizontal')
    ax[1][2].set_ylim([0.0,4.0])
    
    ax[2][0].plot(n/wr.n0, wr.E_symm(n))
    ax[2][0].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[2][0].set_ylabel(r'$E_{bind}$',fontsize=18, rotation='horizontal') 
    
    ax[2][2].plot(n/wr.n0, map(lambda z: wr.np_eq(z)/z, n))
    ax[2][2].plot(n/wr.n0, array([0.14 for i in n]), color='red', linewidth=2.5)
    ax[2][2].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[2][2].set_ylabel(r'$\frac{n_p(n)}{n}$',fontsize=18, rotation='horizontal') 
    
    ax[2][1].plot(n/wr.n0, map(lambda z: wr.C.eta_o(wr.f_eq(z - wr.np_eq(z),wr.np_eq(z))),n))
    ax[2][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
    ax[2][1].set_ylabel(r'$\eta_\omega(f(n))$',fontsize=18, rotation='horizontal')
    ax[2][1].set_ylim([0.0,4.0])
    
    
        
    
#         ax[0][0].plot(x, map(lambda z: wr.C.phi_n(z),x))
#         ax[0][1].plot(n, wr.P(n))
#         ax[0][2].plot(n, map(lambda z: wr.f_eq(z - wr.np_eq(z),wr.np_eq(z)), n))
#         ax[1][0].plot(x, map(lambda z: wr.C.eta_o(z),x))
#         ax[1][1].semilogy(n_Constr, wr.P_symm(n_Constr))
#         ax[1][2].plot(x, map(lambda z: wr.C.eta_r(z),x)) 
#         ax[2][0].plot(n, wr.E_symm(n)) 
draw()


if masses:
    mmax = []
    #                 
    for i, wr in enumerate(wrlist):
        n_star = arange(wr.n0, 8*wr.n0, wr.n0/10)
        wr.filename='%i.dat'%i
        wr.tabulate()
        M, R = wr.star(wr.n0, 8*wr.n0, wr.n0/10)
        mmax.append(max(M))
        iMax = M.index(max(M))
        nmax = n_star[iMax]
        ax[2][1].plot(n_star/wr.n0, M)
        ax[2][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[2][1].set_ylabel(r'$\frac{M}{M_\odot}$',fontsize=18, rotation='horizontal')
        ax[2][1].annotate('$M_{max} = %.2f$'%max(M), xy=(nmax/wr.n0 ,max(M)), fontsize=10, 
                          arrowprops=dict(facecolor='black',shrink=0.05, width=0.5, frac=0.15), xytext=(nmax, max(M)))
    print mmax


fig.legend(lines, (wr.name for wr in wrlist))
fig.tight_layout()

for wr in wrlist:
    C = wr.C
    print wr.name, C.C_s**2, C.C_o**2, C.C_r**2, C.b, C.c

ioff()
show()




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