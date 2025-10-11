from common_modules import np, sys

from .update import Update
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
            for e_i in e_[1:]: # normalize vectors
                e_i /= np.sqrt(np.dot(e_i, e_i))
            w_ = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])
        else:
            print("Error: D and Q configuration not found", file=sys.stderr)

        self.e_ = e_
        self.w_ = w_  

    # ================= initialize grid quantities ================= #
    def init_physical_properties(self, data: dict):
        keys = ["viscosity", "omega", "c", "u0", "dt", "tau"]  # Define necessary keys
        for key in keys:  
            setattr(self, key, data.get(key))

    def init_grid(self, lx, ly, nx ,ny):
        self.lx = lx  # length in x
        self.ly = ly  # length in y
        self.nx = nx  # number of cells in x
        self.ny = ny  # number of cells in y
        self.wall = np.zeros((ny, nx)) # wall cells in grid

        # define quantities
        self.f = np.zeros((ny, nx, self.Q))

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
        update = Update(self)
        update.stream()
        update.bounce()
        # apply boundary conditions
        if (bc_dict):
            for loc, vals in bc_dict.items():
                self.boundary_condition(loc, vals[0], vals[1])
        update.collide()


        # update.stream()
        # update.bounce()
        # update.collide()

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
 

    # ================= get output quantities ================= # 
    def get_rho(self):
        return self.rho
    
    def get_velocity(self):
        return self.u
    
    def get_kinetic_energy(self):         
        return np.multiply(self.u[:,:,0], self.u[:,:,0]) + np.multiply(self.u[:,:,1], self.u[:,:,1])