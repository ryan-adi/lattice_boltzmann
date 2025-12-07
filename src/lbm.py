from common_modules import np, sys

from .update import *
from .boundary_condition import BoundaryCondition

class LatticeBoltzmann():
    def __init__(self, D, Q):
        # dimension
        configuration = [[1,3], [2,9], [3,27]]
        self.D = D
        self.Q = Q

        _e = np.array([])
        _w = np.array([])
        isInConfiguration = [D,Q] in configuration
        if (isInConfiguration):
            # for D2Q9
            _e = np.array([np.array([0.0, 0.0]),
                            np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                            np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                            np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                            np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
            _w = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])
        else:
            print("Error: D and Q configuration not found", file=sys.stderr)

        # weights and directions
        self._e = _e
        self._w = _w

    # ================= INITIALIZER ================= #
    def init_physical_quantities(self, data: dict):
        keys = ["viscosity", "speedOfSound", "u0", "dt", "density", "relaxationTime"]  # Define necessary keys
        for key in keys:
            setattr(self, key, data.get(key))

    def init_grid_quantities(self, geometry:dict):
        self.lx = geometry["lx"]  # length in x
        self.ly = geometry["ly"]  # length in y
        self.nx = geometry["nx"]  # number of cells in x
        self.ny = geometry["ny"]  # number of cells in y

    def init_field_quantities(self, initFieldConditions:dict):
        nx = self.nx
        ny = self.ny
    
        u0 = initFieldConditions["velocity"]
        rho0 = initFieldConditions["density"]
        random = initFieldConditions["random"]

        # wall
        self.wall = np.zeros((ny, nx),dtype=int) # wall cells in grid

        # Macroscopic density and velocity
        self.rho = rho0 * np.ones((ny, nx))                  # density
        self.vel = np.zeros((ny,nx,self.D))
        self.vel[:,:,0] = u0[0] * np.ones((ny, nx))            # velocity in x
        self.vel[:,:,1] = u0[1] * np.ones((ny, nx))            # velocity in y

        # get equilibrium distribution
        self.f = self.get_f_eq()
        if random:
            self.f += .002 * np.random.randn(ny, nx, self.Q) # induce randomness
    
    # ================= CALC FUNCTIONS ================= #
    def calc_rho(self):
        self.rho = np.sum(self.f, axis=2)

    def calc_vel(self):
        self.vel = calc_fe(self.f, self._e) / self.rho[...,np.newaxis]

    def get_f_eq(self):
        vel = self.vel
        nx = self.nx
        ny = self.ny
        e = self._e
        w = self._w

        f_eq = np.zeros((ny,nx,self.Q))
        for qi, ei, wi in zip(range(self.Q), self._e, self._w):
            cx, cy = ei[0], ei[1]
            temp = wi * self.rho[:,:] * (1 + 3 * (cx*vel[:,:,0] + cy*vel[:,:,1]) 
                                    + (3 * (cx*vel[:,:,0] + cy*vel[:,:,1]))**2 / 2
                                    - 1.5 * (vel[:,:,0]**2+vel[:,:,1]**2))
            f_eq[:,:,qi] = temp 
        return f_eq

    # ================= LBM FUNCTIONS ================= #

    def collide(self):
        # create copy of commonly used variables
        _w = self._w
        _e = self._e
        nx, ny = self.nx, self.ny
        Q = self.Q
        dt = self.dt
        tau = self.relaxationTime

        # update micro velocities
        f_eq = self.get_f_eq()
        self.f +=  dt/tau * (-self.f + f_eq)

    def update(self, bc_dict):
        # stream & bounce
        stream(self.f, self._e)
        bounce(self.f, self.wall)

        # apply boundary conditions
        bc_loc = ["left", "right", "top", "bottom"]
        for loc in bc_loc:
            vals = bc_dict[loc]
            if (vals):
                self.boundary_condition(loc, vals[0], vals[1])

        # apply corner boundary conditions
        if (bc_dict["corner"]):
            self.corner_boundary_condition() 
            
        self.calc_rho()
        self.calc_vel()
        self.collide()

    # ================= boundary conditions ================= #
    def boundary_condition(self, location, type, val):
        bc = BoundaryCondition(self)

        if (type=='velocity'):
            if (location=='left'):
                bc.velocity_left_bc(val)
            elif (location=='right'):
                bc.velocity_right_bc(val)
            elif (location=='top'):
                bc.velocity_top_bc(val)
            elif (location=='bottom'):
                bc.velocity_bottom_bc(val)
            else :
                print("BC location not found", file=sys.stderr)
        elif (type=='pressure'):
            if (location=='left'):
                bc.pressure_left_bc(val)
            # elif (location=='right'):
            #     bc.pressureright_bc(val)
            # elif (location=='top'):
            #     bc.pressuretop_bc(val)
            # elif (location=='bottom'):
            #     bc.pressurebottom_bc(val)
            else :
                print("BC location not found", file=sys.stderr)
        elif (type=='set'): # set micro velocity BC
            if (location=='left'):
                bc.set_left_bc(val)
            elif (location=='right'):
                bc.set_right_bc(val)
            elif (location=='top'):
                bc.set_top_bc(val)
            elif (location=='bottom'):
                bc.set_bottom_bc(val)
            else :
                print("BC location not found", file=sys.stderr)
        elif (type=='absorb'): # set absorbing BC
            if (location=='left'):
                bc.absorb_left_bc()
            elif (location=='right'):
                bc.absorb_right_bc()
            elif (location=='top'):
                bc.absorb_top_bc()
            elif (location=='bottom'):
                bc.absorb_bottom_bc()
            else :
                print("BC location not found", file=sys.stderr)
        else:
            pass
            #print("No boundary condition defined on "+location)

    def corner_boundary_condition(self):
        bc = BoundaryCondition(self)
        bc.corner_bc()

    # ================= get macro quantities ================= #
    def get_rho(self):
        return self.rho

    def get_velocity(self):
        return self.vel

    def get_kinetic_energy(self):
        return np.multiply(self.vel[:,:,0], self.vel[:,:,0]) + np.multiply(self.vel[:,:,1], self.vel[:,:,1])

    def get_curl(self):
        dx_uy = (self.vel[1:-1,2:,1]-self.vel[1:-1,0:-2,1]) / 2
        dy_ux = (self.vel[2:,1:-1,0]-self.vel[0:-2,1:-1,0]) / 2
        return dx_uy - dy_ux