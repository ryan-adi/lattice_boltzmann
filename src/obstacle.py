from common_modules import np

from .lbm import LatticeBoltzmann

class Obstacle():
    def __init__(self, lbm:LatticeBoltzmann):
        self.lbm = lbm

    def create_box(self, bb_min: list[int], bb_max: list[int]):
        '''
        Initialize obstacle cells depending on bounding box definition
        @param bb_min := min bounding box, index of grid cells that are obstacle
        @param bb_max := max bounding box, index of grid cells that are obstacle
        '''
        for yi in range(bb_min[1], bb_max[1]+1):
            for xi in range(bb_min[0], bb_max[0]+1):
                isInsideX = (xi < self.lbm.nx) & (xi > -1)
                isInsideY = (yi < self.lbm.ny) & (yi > -1)
                if (isInsideX & isInsideY):
                    self.lbm.wall[yi, xi] = 1

    def create_circle(self, radius:float, center:np.array):
        '''
        Initialize obstacle cells depending on bounding box definition
        @param radius := radius of circle
        @param center := center of circle
        '''
        ny = self.lbm.ny
        nx = self.lbm.nx
        distance = lambda x, y : np.dot((y-x),(y-x)) 
        for yi in range(ny):
            for xi in range(nx):
                position = np.array([xi,yi])
                if ((distance(center, position)) < radius*radius) : 
                    self.lbm.wall[yi, xi] = 1 

    def move(self, displacement:np.array):
        '''
        Move obstacles in a specific direction
        @param displacement := displacement vector
        '''

        self.lbm.wall[:,:] = np.roll(self.lbm.wall[:,:], displacement[0], axis=1)
        self.lbm.wall[:,:] = np.roll(self.lbm.wall[:,:], displacement[1], axis=0)

        q = self.lbm.Q
        e_ = self.lbm.e_
        # streaming
        for i, ei in zip(range(q), e_):
            cx, cy = ei[0], ei[1]
            self.lbm.f[:,:,i] = np.roll(self.lbm.f[:,:,i], cx, axis=1)
            self.lbm.f[:,:,i] = np.roll(self.lbm.f[:,:,i], cy, axis=0)