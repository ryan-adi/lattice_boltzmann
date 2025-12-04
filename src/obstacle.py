from common_modules import np

from .lbm import LatticeBoltzmann

class Obstacle():
    def __init__(self, lbm:LatticeBoltzmann):
        self.lbm = lbm

    def create_obstacle(self, obstacle_dict):
        shape = obstacle_dict["shape"]
        match shape:
            case "box":
                bb_min = obstacle_dict["bb_min"]
                bb_max = obstacle_dict["bb_max"]
                self.create_box(bb_min, bb_max)
            case "circle":
                radius = obstacle_dict["radius"]
                center = obstacle_dict["center"]
                self.create_circle(radius, center)

    def create_box(self, bb_min: list[int], bb_max: list[int]):
        '''
        Initialize obstacle cells depending on bounding box definition
        @param bb_min := min bounding box, index of grid cells that are obstacle
        @param bb_max := max bounding box, index of grid cells that are obstacle
        '''
        y_min = max(bb_min[1], 0)
        y_max = min(bb_max[1], self.lbm.ny-2)+1
        x_min = max(bb_min[0], 0)
        x_max = min(bb_max[0], self.lbm.nx-2)+1
        # inner cells
        self.lbm.wall[y_min:y_max, x_min:x_max] = 1
        # boundary cells
        self.lbm.wall[y_min-1,x_min-1:x_max+1] = 2
        self.lbm.wall[y_max,x_min-1:x_max+1] = 2
        self.lbm.wall[y_min-1:y_max+1,x_min-1] = 2
        self.lbm.wall[y_min-1:y_max+1,x_max] = 2

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
                # boundary
                if ((distance(center, position)) <= (radius+1)**2) : 
                    self.lbm.wall[yi, xi] = 2
                # inner 
                if ((distance(center, position)) <= radius**2) : 
                    self.lbm.wall[yi, xi] = 1 

    def move(self, displacement:np.array):
        '''
        Move obstacles in a specific direction
        @param displacement := displacement vector
        '''
        q = self.lbm.Q
        _e = self.lbm._e

        # move wall
        self.lbm.wall = np.roll(self.lbm.wall[:,:], displacement[1], axis=0)
        self.lbm.wall = np.roll(self.lbm.wall[:,:], displacement[0], axis=1)
        
        # move f
        for i, ei in zip(range(q), _e):
            self.lbm.f[:,:,i] = np.roll(self.lbm.f[:,:,i], displacement[0], axis=1)
            self.lbm.f[:,:,i] = np.roll(self.lbm.f[:,:,i], displacement[1], axis=0)

        # update f on boundaries
        boundary_idx = self.lbm.wall==2
        for i in range(1,5):
            val = max(0, np.dot(displacement, _e[i]))
            self.lbm.f[:,:,i] += 1.0 * val * boundary_idx
        for i in range(5,q):
            val = max(0, np.dot(displacement, _e[i]))
            self.lbm.f[:,:,i] += 1.0 * val * boundary_idx

        ## collide
        # f = self.lbm.f
        # self.rho = np.sum(f, 2)
        # self.vel = np.einsum("ijk,kl->ijl", f, _e) / self.rho[...,np.newaxis]