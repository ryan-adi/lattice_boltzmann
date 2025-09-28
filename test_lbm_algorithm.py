import unittest
from lbm import *
import numpy as np

class TestStringMethods(unittest.TestCase):
    ## for template / exmample only ##
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

class TestStream(unittest.TestCase):
       
        # fluid properties
        viscosity = 0.002                     # viscosity
        omega = 1./(3*viscosity + 0.5)        # relaxation parameter (a function of viscosity)
        u0 = 0.1 * np.array([1.0, 0.0])  # initial in-flow velocity
        c = 1.0                               # lattice velocity
        #rho_base = 1000                            # base density

        ## LBM PARAMS (CURRENTLY ONLY FOR D=2)
        D = 2
        Q = 9
        # for D2Q9
        e_k = np.array([np.array([0.0, 0.0]),
                        np.array([1.0, 0.0]),
                        np.array([0.0, 1.0]),
                        np.array([-1.0, 0.0]),
                        np.array([0.0, -1.0]),
                        np.array([1.0, 1.0]),
                        np.array([-1.0, 1.0]),
                        np.array([-1.0, -1.0]),
                        np.array([1.0, -1.0])])
        for e_ki in e_k[1:]: # noramlize vectors
            e_ki /= np.sqrt(np.dot(e_ki, e_ki))
        w_k = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])


       
        height = 2
        width = 2
        # Equilbrium distributions
        f = np.zeros((height*width, Q))
        
        for i in range(height*width):
            # init equilibrium distribution
            for qi in range(Q):
                f[i, qi] = w_k[qi] * (1 + 3 * np.dot(u0, e_k[qi]) + 4.5 * (np.dot(u0, e_k[qi]))**2 - 1.5 * np.dot(u0, u0)) 

if __name__ == '__main__':
    unittest.main()