import collections
import eosCore
import os
from numpy import *

def dec_map_method(fun):
    def wrap_args(self,x, C):
        if isinstance(x, collections.Iterable):
            map(lambda z: fun(self,z,C), x)
        else:
            fun(self,x,C)
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
            
            
    def tabulate(self):
        if not os.path.isfile(self.filename):
            n = arange(0.001, 4.0, 0.01)
    
    @dec_map_method
    def P(self, n):
        return eosCore.P(n, C)
    
    
    