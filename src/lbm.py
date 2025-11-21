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
        w_ = np.array([])
        isInConfiguration = [D,Q] in configuration
        if (isInConfiguration):
            # for D2Q9
            _e = np.array([np.array([0.0, 0.0]),
                            np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                            np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                            np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                            np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
            # for ei in _e[1:]: # normalize vectors
            #     ei /= np.sqrt(np.dot(ei, ei))
            w_ = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])
        else:
            print("Error: D and Q configuration not found", file=sys.stderr)

        # weights and directions
        self._e = _e
        self.w_ = w_

    
    # ================= update ================= #
    def update(self, bc_dict):

        # # apply boundary conditions
        # if (bc_dict):
        #     for loc, vals in bc_dict.items():
        #         if (vals):
        #             self.boundary_condition(loc, vals[0], vals[1])

        stream(self.f, self._e)
        bounce(self.f, self.wall)

        # apply boundary conditions
        if (bc_dict):
            for loc, vals in bc_dict.items():
                if (vals):
                    self.boundary_condition(loc, vals[0], vals[1])

        # apply corner boundary conditions
        # self.corner_boundary_condition() 

        self.collide()

    def collide(self):
        # create copy of commonly used variables
        w_ = self.w_
        _e = self._e
        u = self.u
        dt = self.dt
        tau = self.tau

        # calc macro quantities
        self.rho = np.sum(self.f, axis=2)
        u = calc_fe(self.f, _e) / self.rho[...,np.newaxis]
        self.u = u
        self.scalar = self.f[:,:,0]

        f_eq = np.zeros(self.f.shape)
        for i, ei, wi in zip(range(9), _e, w_):
            cx, cy = ei[0], ei[1]
            f_eq[:,:,i] = wi * self.rho * (1 + 3 * (cx*u[:,:,0] + cy*u[:,:,1]) 
                                    + (3 * (cx*u[:,:,0] + cy*u[:,:,1]))**2 /2
                                    - 3 * (u[:,:,0]**2+u[:,:,1]**2)/2)
        self.f +=  - dt/tau * (self.f - f_eq)

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
            print("Boundary condition not defined")

    def corner_boundary_condition(self):
        bc = BoundaryCondition(self)

        bc.corner_bc()

        
    # ================= get macro quantities ================= #
    def get_rho(self):
        return self.rho

    def get_velocity(self):
        return self.u

    def get_kinetic_energy(self):
        return np.multiply(self.u[:,:,0], self.u[:,:,0]) + np.multiply(self.u[:,:,1], self.u[:,:,1])

    def get_curl(self):
        dx_uy = (self.u[1:-1,2:,1]-self.u[1:-1,0:-2,1]) / 2
        dy_ux = (self.u[2:,1:-1,0]-self.u[0:-2,1:-1,0]) / 2
        return dx_uy - dy_ux
    
    def get_scalar(self):
        return self.scalar