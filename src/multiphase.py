from common_modules import np, nb, plt
from src.lbm import *
from src.update import *

class Multiphase(LatticeBoltzmann):
    def __init__(self, D, Q):
        super().__init__(D, Q) 

    # ================= INITIALIZER ================= #
    def init_physical_quantities(self, data: dict):
        keys = ["viscosity", "speedOfSound", "u0", "dt", "density", "relaxationTime"]  # Define necessary keys
        for key in keys:
            setattr(self, key, data.get(key))

    def init_multiphase_components(self, data:dict):
        keys = ["config", "n_fluid", "radius"]
        for key in keys:
            setattr(self, key, data.get(key))

        self.tau = data["relaxationTime"]
        self.g_coeff = np.zeros((self.n_fluid, self.n_fluid))
        for gi, g_coeff in enumerate(data["g_coeff"]):
            x, y = gi//self.n_fluid, gi%self.n_fluid
            self.g_coeff[x,y] = g_coeff

    def init_field_quantities(self, initFieldConditions):
        nx = self.nx
        ny = self.ny
        n_fluid = self.n_fluid
    
        u0=initFieldConditions["velocity"]
        rho0=initFieldConditions["density"]
        random=initFieldConditions["random"]

        # wall
        self.wall = np.zeros((ny, nx),dtype=int) # wall cells in grid

        self.rho = np.zeros((self.n_fluid, self.ny, self.nx))
        self.vel = np.zeros((self.ny, self.nx,2)) 
        match self.config:
            case 1:
                circle = np.zeros((ny,nx))
                for i in range(nx):
                    for j in range(ny):
                        if ((i-nx//2)**2+(j-ny//2)**2 < self.radius**2):
                            circle[j,i] = 1.0
                            
                self.rho[0,:,:] = 2.0 * circle + 0.1
                self.vel[:,:,0] = u0[0] * circle 
                self.vel[:,:,1] = u0[1] * circle 
            case 2:
                # two phase component
                self.rho[0,:ny//2,:] = 0.1
                self.rho[1,:ny//2,:] = 0.9
                self.rho[0,ny//2:,:] = 0.9
                self.rho[1,ny//2:,:] = 0.1
            case 3:
                # equal density start
                self.rho = rho0 * np.ones((n_fluid, ny, nx))                  # density
                self.vel = np.zeros((ny,nx,self.D))
                self.vel[:,:,0] = u0[0] * np.ones((ny, nx))            # velocity in x
                self.vel[:,:,1] = u0[1] * np.ones((ny, nx))            # velocity in y 
                
        self.f = self.get_f_eq()

        self.force = np.zeros((n_fluid,ny,nx,2))               # force

        # get equilibrium distribution
        self.f = self.get_f_eq()
        if random:
            self.f += .002 * np.random.randn(n_fluid, ny, nx, self.Q) # induce randomness

    # ============ GET FUNCTION ============ #
    def get_rho(self):
        return np.sum(self.rho, axis=0) # (ny, nx)
    
    def get_psi(self, component_rho):
        return self.density * (1. - np.exp(-component_rho/self.density)) # (ny, nx)  

    def get_f_eq(self):
        vel = self.vel
        e = self._e
        w = self._w

        f_eq = np.zeros((self.n_fluid, self.ny, self.nx,self.Q))
        for fluid_i in range(self.n_fluid):
            for qi, ei, wi in zip(range(self.Q), self._e, self._w):
                cx, cy = ei[0], ei[1]
                f_eq[fluid_i,:,:,qi] = wi * self.rho[fluid_i,:,:] * (1 + 3 * (cx*vel[:,:,0] + cy*vel[:,:,1]) 
                                        + (3 * (cx*vel[:,:,0] + cy*vel[:,:,1]))**2 / 2
                                        - 1.5 * (vel[:,:,0]**2+vel[:,:,1]**2))
                
        return f_eq   
    
    # ============ CALC FUNCTIONS ============ #
    def calc_rho(self):
        self.rho = np.sum(self.f, axis=3)

    def calc_vel(self):
        self.vel[:,:,:] = 0.
        for fluid_i in range(self.n_fluid):
            vel = calc_fe(self.f[fluid_i,:,:,:], self._e) #/ self.rho[...,np.newaxis]
            self.vel[:,:,0] += vel[:,:,0] + self.force[fluid_i,:,:,0] * 0.5
            self.vel[:,:,1] += vel[:,:,1] + self.force[fluid_i,:,:,1] * 0.5
        
        rhoTotal = self.get_rho()
        self.vel[:,:,0] /= rhoTotal
        self.vel[:,:,1] /= rhoTotal

    def calc_sc_forces(self):
        nx = self.nx
        ny = self.ny
        n_fluid = self.n_fluid

        for fluid_i in range(n_fluid):
            self.force[fluid_i,:,:,0] = 0.0
            self.force[fluid_i,:,:,1] = 0.0
            for fluid_j in range(n_fluid):
                fx_temp = np.zeros((ny,nx))
                fy_temp = np.zeros((ny,nx))
                for ei, wi in zip(self._e[1:,:], self._w[1:]):
                    cx, cy = int(ei[0]), int(ei[1])

                    temp = self.rho[fluid_j,:,:].copy()                   
                    cx, cy = int(ei[0]), int(ei[1])
                    rho_shift = np.roll(temp[:,:], -cx, axis=1)
                    rho_shift = np.roll(rho_shift, -cy, axis=0)
                    
                    psinb = self.get_psi(rho_shift)
                    fx_temp[:,:] += wi * ei[0] * psinb
                    fy_temp[:,:] += wi * ei[1] * psinb
                
                psiloc = self.get_psi(self.rho[fluid_i,:,:])
                fx_temp[:,:] *= - self.g_coeff[fluid_i,fluid_j] * psiloc
                fy_temp[:,:] *= - self.g_coeff[fluid_i,fluid_j] * psiloc

                self.force[fluid_i,:,:,0] += fx_temp
                self.force[fluid_i,:,:,1] += fy_temp

    # ============ LBM FUNCTIONS ============ #
    def collide(self):
        vel = self.vel
        force = self.force

        # compute equilibrium
        f_eq = self.get_f_eq()
        
        for fluid_i in range(self.n_fluid):
            omega = self.dt / self.tau[fluid_i]

            forcing = np.zeros(((self.ny, self.nx, self.Q)))
            # compute forcing
            for qi, ei, wi in zip(range(self.Q), self._e, self._w):
                cx, cy = ei[0], ei[1]
                x_term =  (3. * (cx - vel[:,:,0]) + 9. * cx * (cx*vel[:,:,0] + cy*vel[:,:,1])) * force[fluid_i,:,:,0]
                y_term =  (3. * (cy - vel[:,:,1]) + 9. * cy * (cx*vel[:,:,0] + cy*vel[:,:,1])) * force[fluid_i,:,:,1]
                forcing[:,:,qi] = wi * (1-0.5*omega) * (x_term + y_term)

            # update micro velocities
            self.f[fluid_i, :, :, :] +=  omega * (-self.f[fluid_i,:,:,:] + f_eq[fluid_i,:,:,:]) + forcing

    def stream_mp(self):
        for fluid_i in range(self.n_fluid):
            stream(self.f[fluid_i,:,:,:], self._e)

    def bounce_mp(self):
        for fluid_i in range(self.n_fluid):
            bounce(self.f[fluid_i,:,:,:], self.wall)

    def boundary_mp(self, bc_dict):
        # apply boundary conditions
        if (bc_dict):
            for loc, vals in bc_dict.items():
                if (vals):
                    self.boundary_condition(loc, vals[0], vals[1])

        # apply corner boundary conditions
        # self.corner_boundary_condition() 

    def update(self, bc_dict):
        self.stream_mp()
        self.bounce_mp()
        self.boundary_mp(bc_dict)
       
        self.calc_rho()
        self.calc_sc_forces()
        self.calc_vel()
        self.collide()