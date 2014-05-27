import eosCore as eos
from pylab import *
from RMF.Wrapper import Wrapper
from numpy import* 
from numpy import max
from gtk._gtk import Orientation
from scipy.misc.common import derivative

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

f0 = 0.26

Cs.set_f0(f0)
Cs.set_gp_beta(2.8)
Cs.set_gp_z(0.65)

Cs.set_gp_alpha(1.0)
Cs.set_gp_gamma(2.0)
Cs.set_gp_w(1.0)

Cs.set_gp_e(0)
Cs.set_gp_d(0)
Cs.set_gp_a(0.0)
Cs.set_gp_fcut(0.22)
acs = 1.8

Csj.set_gp_alpha(1.0)
Csj.set_gp_z(0.65)
Csj.set_gp_gamma(1.4)
Csj.set_gp_beta(2.8)
Csj.set_gp_w(1.0)
Csj.set_gp_e(-0.0)

Csj.set_f_tilde(10.0)

#--------------------------------
acj = 0.85
Csj.set_gp_d(-10.0)

Csj.set_omega_a(-12.0)
Csj.set_fcut_omega(0.55)
Csj.set_omega_c(-400.0)

Csj.set_gp_a(3.5)
Csj.set_gp_fcut(0.35)

#---------------------------------

Csj.set_omega_b(0.0)

Csj.set_rho_a(-25.0)
Csj.set_rho_b(0.0)
Csj.set_fcut_rho(0.25)

Csj.set_fcut_omega_l(0.17)
Csj.set_omega_a_l(0*20.0)
Csj.set_fcut_sigma_l(0.17)
Csj.set_sigma_a_l(0*94.88)
Csj.set_fcut_phi_l(0.17)
Csj.set_phi_a_l(0.0)
#------------------------------

acj2 = 0.85
Csj2.set_gp_d(-10.0)

Csj2.set_omega_a(-12.0)
Csj2.set_fcut_omega(0.55)
Csj2.set_omega_c(-400.0)

Csj2.set_gp_a(3.5)
Csj2.set_gp_fcut(0.35)
#----------------------------

Csj2.set_gp_z(0.65)
Csj2.set_gp_gamma(3.0)
Csj2.set_gp_beta(2.8)
Csj2.set_gp_w(1.0)
Csj2.set_gp_e(0.0)
Csj2.set_f_tilde(0.6)
Csj2.set_rho_a(0.0)
Csj2.set_rho_b(10.0)
Csj2.set_omega_b(0.0)
Csj2.set_phi_c(0.0)
Csj2.set_phi_cpow(3.0)
Csj2.set_phi_cfcut(0.65)
Csj2.set_fcut_omega_l(0.1)
Csj2.set_omega_a_l(0.0)
Csj2.set_fcut_sigma_l(0.1)
Csj2.set_sigma_a_l(0.0)
masses = 0


f_tilde_list = [1.5]
alist = [acj, acj2]
wr = Wrapper(C, './tab/rho_test.dat')
wr0 = Wrapper(C0, './tab/rho_test.dat')

wrs = Wrapper(Cs, '')
wrj = Wrapper(Csj,'')
wrj2 = Wrapper(Csj2,'')
wrKV.name = 'KVOR'
wrj.name = 'I'
wrj2.name = 'II'
wrlist = [wrj, wrj2, wrKV]

x = arange(0.0, 1.0, 0.01)
n = arange(0.01, 8*wr.n0, wr.n0/10)



a_seq = arange(1.0, 1.2, 0.05)

for f in linspace(0.195, f0, 5):
    for wr in wrlist:
        if wr != wrKV:
            wr.C.set_f0(f)
            wr.Solve()

for a in linspace(1.0, acs, 3):
    print 'acs=',a
    wrs.C.set_gp_alpha(a)
    wrs.Solve()
  
for a in linspace(1.0, acj,3):
    print 'acj=', a
    wrj.C.set_gp_alpha(a)
    wrj.Solve()
    
for a in linspace(1.0, acj2, 3):
    print 'acj=', a
    wrj2.C.set_gp_alpha(a)
    wrj2.Solve()

mn = 4.78
print (Csj.C_o**2/Csj.eta_o(0) - Csj.C_s**2/Csj.eta_s(0))/4.78**2
print (KV.C_o**2/KV.eta_o(0) - KV.C_s**2/KV.eta_s(0))/4.78**2
print wrj.Asymm(wrj.n0)
def testAPR():
    cAPR = 1.0/0.16
    n = arange(0.001, 3.0, 0.01)
    nAPR=cAPR*array([0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
    APR = array([-6.48, -12.13, -15.04, -16.00, -15.09, -12.88, -5.03, 2.13, 15.46, 34.39, 58.35, 121.25, 204.02])
    nAPR_N = cAPR*array([0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.8, 0.96])
    APR_N = array([4.45, 6.45, 9.65, 13.29, 17.94, 22.92, 27.49, 38.82, 54.95, 75.13, 99.75, 127.58, 205.34, 305.87])
    
    nGN = cAPR*array([0.08,0.12,0.16,0.20,0.24, 0.28,  0.32, 0.36,0.40,0.44,0.48, 0.52,  0.56    ])
    GN = array([ 8.76, 11.93,16.06,20.93,27.02,33.83,41.13,49.65,58.65,68.39,79.59,91.03,103.69]) 
    
    nG = cAPR*array([0.08,0.10,0.12, 0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32, 0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58])
    G = array([-11.93, -13.87, -14.85, -15.82, -16.31, -15.82, -15.33, -14.36, -13.14, -11.93, -10.2, -8.52, -6.57, -4.38, -1.7, 0.49, 3.41, 6.57, 9.49, 12.9, 16.31, 20.20, 24.34, 28.23, 32.86, 37.48 ])
    #lapr, = plot(nAPR, APR)
    lG, = plot(nG, G)
    #plot(nAPR_N, APR_N, color=lapr.get_color())
    plot(nGN, GN, color=lG.get_color())
    xlabel(r'$\frac{n}{n_0}$',fontsize=18)
    ylabel(r'$E(n)$', rotation='horizontal', fontsize = 18)
    #l1, = plot(n/wrj.n0, wrj.E_symm(n))
    #plot(n/wrj.n0, wrj.E_n(n), color=l1.get_color())
    l2, = plot(n/wrj.n0, wrKV.E_symm(n))
    plot(n/wrj.n0, wrKV.E_n(n), color=l2.get_color())
    #legend([lapr, lG, l1, l2],['APR','Gandolfi', 'II','KVOR'],loc=2)
    legend([lG, l2],['Gandolfi', 'II'],loc=2)
    show()

def testAsymm():
    nAsymm = arange(0.01, 0.75, 0.01)
    c = 1.0/0.16
    nLow = c*array([0.05, 0.1, 0.15])
    nUp = c*array([0.05, 0.1, 0.15])
    yLow = array([12.0, 22.0, 26.0])
    yUp = array([16.0, 24.0, 34.0])
    plot(nLow, yLow, nUp, yUp, c='r', linewidth = 2.0)
    l1,= plot(nAsymm/wrKV.n0, wrKV.Asymm(nAsymm))
    l2,=plot(nAsymm/wrj.n0, wrj.Asymm(nAsymm))
    xlabel(r'$\frac{n}{n_0}$',fontsize=18)
    ylabel(r'$A_{symm}(n)$', rotation='horizontal', fontsize = 18)
    legend([ l1, l2],['KVOR','II'],loc=2)
    show()
    
# testAsymm()
    
testAPR()
    
def Vs(wr, n):
    return derivative(wr.P, n, dx=1e-1,order = 5)/derivative(wr.E, n, dx=1e-1,order = 5)

# plot(n, wrKV.P(n))
# show()

# plot(arange(0.3, 4.0, 0.1), map(lambda z:Vs(wrKV,z),arange(0.3, 4.0, 0.1)))
# show()

#testAPR() 

ion()
fig, ax = subplots(3,3)

fig.set_dpi(100)
    
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
def plot1():
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
        
#         ax[2][0].plot(n/wr.n0, wr.E_symm(n))
#         ax[2][0].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
#         ax[2][0].set_ylabel(r'$E_{bind}$',fontsize=18, rotation='horizontal') 
        
        ax[2][0].plot(x, map(lambda z: wr.C.eta_r(z),x))
        ax[2][0].set_xlabel(r'$f$',fontsize=18, rotation='horizontal')
        ax[2][0].set_ylabel(r'$\eta_\sigma(f)$',fontsize=18, rotation='horizontal')
        ax[2][0].set_ylim([0.0,4.0])
        
        ax[2][2].plot(n/wr.n0, map(lambda z: wr.np_eq(z)/z, n))
        ax[2][2].plot(n/wr.n0, array([0.14 for i in n]), color='red', linewidth=2.5)
        ax[2][2].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[2][2].set_ylabel(r'$\frac{n_p(n)}{n}$',fontsize=18, rotation='horizontal') 
        


def plot2():
    for wr in wrlist:
        feq = map(lambda z: wr.f_eq(z - wr.np_eq(z), wr.np_eq(z)),n)
        line, = ax[0][0].plot(n/wr.n0, map(lambda z: wr.C.phi_n(z),feq))
        lines.append(line)
        ax[0][0].set_xlabel(r'$n/n0$', fontsize = 18)
        ax[0][0].set_ylabel(r'$\Phi(n)$', fontsize=18, rotation='horizontal' )
        
        ax[0][1].plot(n/wr.n0, wr.P(n))
        ax[0][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[0][1].set_ylabel(r'$P_{NSM}(n)$',fontsize=18, rotation='horizontal')
        
        ax[0][2].plot(n/wr.n0, map(lambda z: wr.C.phi_n(wr.f_eq(z - wr.np_eq(z),wr.np_eq(z))), n))
        ax[0][2].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[0][2].set_ylabel(r'$\frac{m^\ast}{m_n}$',fontsize=18, rotation='horizontal')
        
        ax[1][0].plot(n/wr.n0, map(lambda z: wr.C.eta_o(z),feq))
        ax[1][0].set_xlabel(r'$n/n0$',fontsize=18, rotation='horizontal')
        ax[1][0].set_ylabel(r'$\eta_\omega(n)$',fontsize=18, rotation='horizontal') 
        
        ax[1][1].semilogy(n_Constr/wr.n0, wr.P_symm(n_Constr))
        ax[1][1].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[1][1].set_ylabel(r'$P_{SNM}(n)$',fontsize=18, rotation='horizontal')
        ax[1][1].set_ylim([1.0, 400.0])
        
        ax[1][2].plot(n/wr.n0, map(lambda z: wr.C.eta_r(z),feq))
        ax[1][2].set_xlabel(r'$n$',fontsize=18, rotation='horizontal')
        ax[1][2].set_ylabel(r'$\eta_\rho(n)$',fontsize=18, rotation='horizontal')
        ax[1][2].set_ylim([0.0,4.0])
        
#         ax[2][0].plot(n/wr.n0, wr.E_symm(n))
#         ax[2][0].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
#         ax[2][0].set_ylabel(r'$E_{bind}$',fontsize=18, rotation='horizontal') 
#         
        ax[2][0].plot(n/wr.n0, map(lambda z: wr.C.eta_s(z),feq))
        ax[2][0].set_xlabel(r'$n$',fontsize=18, rotation='horizontal')
        ax[2][0].set_ylabel(r'$\eta_\sigma(n)$',fontsize=18, rotation='horizontal')
        ax[2][0].set_ylim([0.0,4.0])

        ax[2][2].plot(n/wr.n0, map(lambda z: wr.np_eq(z)/z, n))
        ax[2][2].plot(n/wr.n0, array([0.14 for i in n]), color='red', linewidth=2.5)
        ax[2][2].set_xlabel(r'$\frac{n}{n_0}$',fontsize=18, rotation='horizontal')
        ax[2][2].set_ylabel(r'$\frac{n_p(n)}{n}$',fontsize=18, rotation='horizontal')

    
plot1()
    
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
    print wr.name, C.C_s, C.C_o, C.C_r, C.b, C.c

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