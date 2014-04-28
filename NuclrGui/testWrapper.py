
import unittest
import sys
import os
print 'ldpath=',os.environ.get('LD_LIBRARY_pPATH')
print ':'.join(x for x in sys.path if x)
from RMF.Wrapper import Wrapper
import eosCore
from numpy import *

C = eosCore.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
wr = Wrapper(C, 'test.dat')
n = arange(0.001, 4.0, 0.01)
class Test(unittest.TestCase):

    
    def testName(self):
        pass
    
    def testP(self):
        print wr.P(n)
        
    def testE(self):
        print wr.E(n)
        
    def testNpeq(self):
        print wr.np_eq(n)
        
    def testFeq(self):
        np_eq = wr.np_eq(n)
        feq = []
        for i, item in enumerate(n):
            feq.append(wr.f_eq(item - np_eq[i], np_eq[i]))
        print feq
    
    def testEsymm(self):
        print wr.E_symm(wr.n0)

    def testEn(self):
        print wr.E_n(wr.n0)
        
    def testPsymm(self):
        print wr.P_symm(wr.n0)
        
    def testPn(self):
        print wr.P_n(wr.n0)
        
    def testSolve(self):
        print wr.C.f0
        wr.Solve()
        
    def testTab(self):
        wr.tabulate()
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()