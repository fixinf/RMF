
import unittest
from RMF import Wrapper
import eosCore
from numpy import *

class Test(unittest.TestCase):
    
    def __init__(self):
        super(Test).__init__()
        self.C = eosCore.KVOR(sqrt(179.56),sqrt(87.6),sqrt(100.64),7.7346e-3,3.4462e-4,0.65)
        self.wr = Wrapper(self.C, 'test.dat')
    
    def testName(self):
        pass
    
    def testP(self):
        n = arange(0.001, 4.0, 0.01)
        print self.wr.P(n, self.C)
        
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()