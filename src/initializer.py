from common_modules import np

from .lbm import LatticeBoltzmann

class Initializer():
    def __init__(self, lbm:LatticeBoltzmann):
        self.lbm = lbm

    def physical_quantities(self, data: dict):
        keys = ["viscosity", "c", "u0", "rho", "dt", "tau"]  # Define necessary keys
        for key in keys:
            setattr(self.lbm, key, data.get(key))

    def grid_quantities(self, lx, ly, nx ,ny):
        self.lbm.lx = lx  # length in x
        self.lbm.ly = ly  # length in y
        self.lbm.nx = nx  # number of cells in x
        self.lbm.ny = ny  # number of cells in y

    def micro_velocities(self, e:np.array, val:np.array):
        nx = self.lbm.nx
        ny = self.lbm.ny

        self.lbm.f = np.ones((ny, nx, self.lbm.Q)) + .01 * np.random.randn(ny, nx, self.lbm.Q)
        for e_i, val_i in zip(e, val):
            self.lbm.f[:,:,e_i] = val_i

    def macro_quantities(self):

        nx = self.lbm.nx
        ny = self.lbm.ny

        self.lbm.wall = np.zeros((ny, nx),dtype=int) # wall cells in grid

        # Macroscopic density and velocity
        self.lbm.rho = np.ones((ny, nx))              # density
        self.lbm.u = np.zeros((ny, nx, self.lbm.D))   # velocity

        # for yi in range(self.ny):
        #     for xi in range(self.nx):

        #         # init equilibrium distribution
        #         for qi in range(self.lbm.Q):
        #             self.f[yi, xi, qi] = self.rho[yi,xi] * self.w_[qi] * (
        #                 1 + 3 * np.dot(self.u0, self.e_[qi]) 
        #                 + 4.5 * (np.dot(self.u0, self.e_[qi]))**2 
        #                 - 1.5 * np.dot(self.u0, self.u0))

        #         # macroscopic values
        #         self.rho[yi, xi] = np.sum(self.f[yi, xi,:])
        #         self.u[yi,xi,:] = np.dot(self.f[yi,xi,1:], self.e_[1:]) * (1-(self.rho[yi,xi]-1)+((self.rho[yi,xi]-1)**2))
    
