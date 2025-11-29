from common_modules import np

from .lbm import LatticeBoltzmann

class Initializer():
    def __init__(self, lbm:LatticeBoltzmann):
        self.lbm = lbm

    def physical_quantities(self, data: dict):
        keys = ["viscosity", "speedOfSound", "u0", "dt", "density", "relaxationTime"]  # Define necessary keys
        for key in keys:
            setattr(self.lbm, key, data.get(key))

    def grid_quantities(self, geometry):
        self.lbm.lx = geometry["lx"]  # length in x
        self.lbm.ly = geometry["ly"]  # length in y
        self.lbm.nx = geometry["nx"]  # number of cells in x
        self.lbm.ny = geometry["ny"]  # number of cells in y

    def field_quantities(self, initFieldConditions):
        nx = self.lbm.nx
        ny = self.lbm.ny

        e=initFieldConditions["microVelDir"]
        val=initFieldConditions["microVelVal"]
        u0=initFieldConditions["velocity"]
        rho0=initFieldConditions["density"]

        self.lbm.f = np.ones((ny, nx, self.lbm.Q)) + .005 * np.random.randn(ny, nx, self.lbm.Q)
        for e_i, val_i in zip(e, val):
            self.lbm.f[:,:,int(e_i)] = val_i

        # wall
        self.lbm.wall = np.zeros((ny, nx),dtype=int) # wall cells in grid

        # Macroscopic density and velocity
        self.lbm.rho = rho0 * np.ones((ny, nx))                  # density
        self.lbm.u = np.zeros((ny,nx,self.lbm.D))
        self.lbm.u[:,:,0] = u0[0] * np.ones((ny, nx))            # velocity in x
        self.lbm.u[:,:,1] = u0[1] * np.ones((ny, nx))            # velocity in y

        # macroscopic scalar
        self.lbm.scalar = np.ones((ny, nx))