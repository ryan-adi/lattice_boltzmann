from common_modules import np

from .lbm import LatticeBoltzmann

class Obstacle():
    def __init__(self, lbm:LatticeBoltzmann):
        self.lbm = lbm

    def create_box(self, bb_min: list[int], bb_max: list[int]):
        '''
        Initialize wall cells depending on bounding box definition
        @param bb_min := min bounding box, index of grid cells that are walls
        @param bb_max := max bounding box, index of grid cells that are walls
        '''
        for yi in range(bb_min[1], bb_max[1]+1):
            for xi in range(bb_min[0], bb_max[0]+1):
                isInsideX = (xi < self.lbm.nx) & (xi > -1)
                isInsideY = (yi < self.lbm.ny) & (yi > -1)
                if (isInsideX & isInsideY):
                    self.lbm.wall[yi, xi] = 1

    def create_circle(self, radius:float, center:np.array):
        '''
        Initialize wall cells depending on bounding box definition
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