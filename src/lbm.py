from common_modules import np, sys

from .update import *
from .boundary_condition import BoundaryCondition

class LatticeBoltzmann():
    def __init__(self, D, Q):
        # initialize
        configuration = [[1,3], [2,9], [3,27]]
        self.D = D
        self.Q = Q

        e_ = np.array([])
        w_ = np.array([])
        isInConfiguration = [D,Q] in configuration
        if (isInConfiguration):
            # for D2Q9
            e_ = np.array([np.array([0.0, 0.0]),
                            np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                            np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                            np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                            np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
            # for e_i in e_[1:]: # normalize vectors
            #     e_i /= np.sqrt(np.dot(e_i, e_i))
            w_ = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])
        else:
            print("Error: D and Q configuration not found", file=sys.stderr)

        self.e_ = e_
        self.w_ = w_

    # ================= initialize grid quantities ================= #
    def init_quantities(self, data: dict):
        keys = ["viscosity", "omega", "c", "u0", "rho", "dt", "tau"]  # Define necessary keys
        for key in keys:
            setattr(self, key, data.get(key))

    def init_grid(self, lx, ly, nx ,ny):
        self.lx = lx  # length in x
        self.ly = ly  # length in y
        self.nx = nx  # number of cells in x
        self.ny = ny  # number of cells in y
        self.wall = np.zeros((ny, nx)) # wall cells in grid

        # define quantities
        self.f = np.ones((ny, nx, self.Q))

        # Macroscopic density and velocity
        self.rho = np.ones((ny, nx))    # density
        self.u = np.zeros((ny, nx, self.D))   # velocity

        for yi in range(self.ny):
            for xi in range(self.nx):

                # init equilibrium distribution
                for qi in range(self.Q):
                    self.f[yi, xi, qi] = self.w_[qi] * (
                        1 + 3 * np.dot(self.u0, self.e_[qi]) + 4.5 * (np.dot(self.u0, self.e_[qi]))**2 - 1.5 * np.dot(self.u0, self.u0))

                # macroscopic values
                self.rho[yi, xi] = np.sum(self.f[yi, xi,:])
                self.u[yi,xi,:] = np.dot(self.f[yi,xi,1:], self.e_[1:]) * (1-(self.rho[yi,xi]-1)+((self.rho[yi,xi]-1)**2))

    def init_wall(self, bb_min: list[int], bb_max: list[int]):
        '''
        Initialize wall cells depending on bounding box definition
        @param bb_min := min bounding box, index of grid cells that are walls
        @param bb_max := max bounding box, index of grid cells that are walls
        '''
        for yi in range(bb_min[1], bb_max[1]+1):
            for xi in range(bb_min[0], bb_max[0]+1):
                isInsideX = (xi < self.nx-1) & (xi > -1)
                isInsideY = (yi < self.ny-1) & (yi > -1)
                if (isInsideX & isInsideY):
                    self.wall[yi, xi] = 1


    # ================= update ================= #
    def update(self, bc_dict):
        self.f = stream(self.f)
        self.f = bounce(self.f, self.wall)

        # apply boundary conditions
        if (bc_dict):
            for loc, vals in bc_dict.items():
                if (vals):
                    self.boundary_condition(loc, vals[0], vals[1])

        self.collide()


    def collide(self):
        # create copy of commonly used variables
        omega = self.omega
        w_ = self.w_
        e_ = self.e_
        u = self.u
        dt = self.dt
        tau = self.tau

        # make cxs and cys
        cxs = e_[:,0]
        cys = e_[:,1]

        # store original f field
        f = self.f

        if (False):
            # calc macro quantities
            rho = np.sum(f, 2)
            ux = np.sum(f*cxs, 2) / rho
            uy = np.sum(f*cys, 2) / rho

            f_eq = np.zeros(f.shape)
            for i, ei, wi in zip(range(9), e_, w_):
                cx, cy = ei[0], ei[1]
                # ux = self.u[:,:,0]
                # uy = self.u[:,:,1]
                f_eq[:,:,i] = self.rho * wi * (1 + 3 * (cx*ux + cy*uy) 
                                        + (3 * (cx*ux + cy*uy))**2 /2
                                        - 3 * (ux**2+uy**2)/2)
            self.f +=  - dt/tau * (self.f - f_eq)

            # store in self
            self.rho = rho
            self.u[:,:,0] = ux
            self.u[:,:,1] = uy

        if (True):
            n_offset = 1
            for yi in range(n_offset, self.ny-n_offset):
                for xi in range(n_offset, self.nx-n_offset):

                    # Skip over cells containing barriers
                    if (self.wall[yi, xi]):
                        continue
                    else:
                    
                        # Compute the macroscopic density
                        self.rho[yi, xi] = np.sum(self.f[yi,xi,:])

                        # Compute the macroscopic velocities
                        if (self.rho[yi,xi] > 0):
                            self.u[yi,xi,:] = np.dot(self.f[yi,xi,1:], e_[1:]) * (1-(self.rho[yi,xi]-1)+((self.rho[yi,xi]-1)**2))

                        # Compute f
                        for qi in range(1,9):
                            # corresponds to cs=0.577
                            f_eq = w_[qi] * (1 + 3 * np.dot(u[yi,xi,:], e_[qi]) + 4.5 * (np.dot(u[yi,xi,:], e_[qi]))**2
                                             - 1.5 * np.dot(u[yi,xi,:], u[yi,xi,:]))
                            self.f[yi, xi, qi] =  (1-dt/tau) * f[yi, xi, qi] + (dt/tau) * f_eq
                        self.f[yi, xi, 0]  = self.rho[yi, xi] - np.sum(f[yi, xi, 1:])

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
            #     bc.pressure_right_bc(val)
            # elif (location=='top'):
            #     bc.pressure_top_bc(val)
            # elif (location=='bottom'):
            #     bc.pressure_bottom_bc(val)
            else :
                print("BC location not found", file=sys.stderr)
        elif (type=='bounce'):
            if (location=='left'):
                bc.bounce_back_left_bc()
            elif (location=='right'):
                bc.bounce_back_right_bc()
            elif (location=='top'):
                bc.bounce_back_top_bc()
            elif (location=='bottom'):
                bc.bounce_back_bottom_bc()
            else :
                print("BC location not found", file=sys.stderr)
        else: # set velocity BC
            if (location=='left'):
                bc.set_velocity_left_bc(val)
            elif (location=='right'):
                bc.set_velocity_right_bc(val)
            elif (location=='top'):
                bc.set_velocity_top_bc(val)
            elif (location=='bottom'):
                bc.set_velocity_bottom_bc(val)
            else :
                print("BC location not found", file=sys.stderr)

        ## for corners
        # bc.no_slip()
        


    # ================= get macro quantities ================= #
    def get_rho(self):
        return self.rho

    def get_velocity(self):
        return self.u

    def get_kinetic_energy(self):
        return np.multiply(self.u[:,:,0], self.u[:,:,0]) + np.multiply(self.u[:,:,1], self.u[:,:,1])
