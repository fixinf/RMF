import collections
import os
from numpy import *
import eosCore


def dec_map_method(fun):
    def wrap_args(self, x):
        if isinstance(x, collections.Iterable):
            return map(lambda z: fun(self, z), x)
        else:
            return fun(self, x)
    return wrap_args

class Wrapper:
    '''
    RMF eosCore functions wrapper
    '''


    def __init__(self, C, filename, needsTab = False):
        '''
        Constructor
        '''
        self.C = C
        self.filename = filename
        if needsTab:
            self.tabulate
            
        self.m_pi = 134.92
        self.Dim = self.m_pi**3
        self.m_n = 938.0
        self.n0 = (197.33/self.m_pi)**3 * 0.16


            
            
    def tabulate(self):
        self.Ctab = eosCore.EoS_tab(self.C, self.filename, 0.001, 8*self.n0, 0.01*self.n0)
    
    @dec_map_method
    def P(self, n):
        return eosCore.P(n, self.C)
    
    @dec_map_method
    def E(self, n):
        return eosCore.E(n, self.C)
    
    @dec_map_method
    def np_eq(self, n):
        return eosCore.np_eq(n, self.C)
    
    def f_eq(self, nn, np):
        return eosCore.f_eq(nn, np, self.C)
    
    def star(self, start, stop, step):
        R = []
        M = []
        for x in arange(start, stop, step):
            r, m, mgr = eosCore.star(x, self.Ctab)
            R.append(r)
            M.append(m)
        return R, M
    
    @dec_map_method
    def E_symm(self, n):
        return eosCore.t_E(n/2, n/2, self.C)/(self.Dim*n) - self.m_n
    
    @dec_map_method
    def E_n(self, n):
        return eosCore.t_E(n, 0, self.C)/(self.Dim*n) - self.m_n
    
    @dec_map_method
    def P_symm(self, n):
        const = 197.33*(0.16/self.n0)**(4.0/3.0)*self.m_pi**(-4.0)
        fun = lambda z: eosCore.t_E(z/2, z/2, self.C)
        dn = 1e-6
        dE = (fun(n + dn) - fun(n))/dn
        return (n*dE - fun(n))*const
    
    @dec_map_method
    def P_n(self, n):
        const = 197.33*(0.16/self.n0)**(4.0/3.0)*self.m_pi**(-4.0)
        fun = lambda z: eosCore.t_E(z, 0.0, self.C)
        dn = 1e-6
        dE = (fun(n + dn) - fun(n))/dn
        return (n*dE - fun(n))*const 
    
    def Solve(self, E0=-16.0, K0=275.0, A0=32.0):
        eosCore.solve(self.C, self.n0, self.C.f0, E0, K0, A0)
    
    
    