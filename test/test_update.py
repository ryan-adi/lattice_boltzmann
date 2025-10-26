from common_modules import unittest, np

from src.update import *

class TestUpdateFunctions(unittest.TestCase):
    def setUp(self):
        # init directions
        self.e_ = np.array([np.array([0.0, 0.0]),
                            np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                            np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                            np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                            np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
        
        # init micro velocities
        self.ny = self.nx = 5
        self.q = 9
        self.f = np.zeros((self.ny,self.nx,self.q))
        for j in range(self.ny):
            for i in range(self.nx):
                self.f[j,i,:] = i+j

        # init walls
        self.wall = np.zeros((self.ny,self.nx))
        self.wall[2,2] = 1

    def test_stream_size(self):
        stream(self.f, self.e_)
        self.assertEqual(self.f.shape, (self.ny,self.nx,self.q))

    def test_stream(self):
        stream(self.f, self.e_)
        self.assertAlmostEqual(self.f[2,2,0], 4)
        self.assertAlmostEqual(self.f[2,2,1], 3)
        self.assertAlmostEqual(self.f[2,2,2], 3)
        self.assertAlmostEqual(self.f[2,2,3], 5)
        self.assertAlmostEqual(self.f[2,2,4], 5)
        self.assertAlmostEqual(self.f[2,2,5], 2)
        self.assertAlmostEqual(self.f[2,2,8], 4)

    def test_bounce_size(self):
        bounce(self.f, self.wall)
        self.assertEqual(self.f.shape, (self.ny,self.nx,self.q))

    def test_bounce(self):
        bounce(self.f, self.wall)
        self.assertAlmostEqual(self.f[2,2,0], 0)
        self.assertNotAlmostEqual(self.f[2,3,0], 0)
        self.assertAlmostEqual(self.f[3,2,2], 4)


if __name__ == '__main__':
    unittest.main()