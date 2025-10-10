import numpy as np

class Update():
    def __init__(self, lb):
        self.lb = lb

    def stream(self):
        f = self.lb.f
        
        # Stream all internal cells
        n_offset = 0 # TODO
        for yi in range(n_offset, self.lb.ny-n_offset):
            for xi in range(n_offset, self.lb.nx-n_offset):

                self.lb.f[yi, xi, 1] = f[yi, xi-1, 1]
                self.lb.f[yi, xi, 2] = f[yi-1, xi, 2]
                self.lb.f[yi, xi-1, 3] = f[yi, xi, 3]
                self.lb.f[yi-1, xi, 4] = f[yi, xi, 4]
                self.lb.f[yi, xi, 5] = f[yi-1, xi-1, 5]
                self.lb.f[yi, xi-1, 6] = f[yi-1, xi, 6]
                self.lb.f[yi-1, xi-1, 7] = f[yi, xi, 7]
                self.lb.f[yi-1, xi, 8] = f[yi, xi-1, 8]
                
        # Tidy up the edges (TODO)
        # xi = self.lb.nx
        # for yi in range(1, self.lb.ny-1):
        #     self.lb.f[yi, xi, 2] = f[yi-1, xi, 2]
        #     self.lb.f[yi-1, xi, 4] = f[yi, xi, 4]


    def collide(self):
        # create copy of commonly used variables 
        omega = self.lb.omega
        w_ = self.lb.w_
        e_ = self.lb.e_
        u = self.lb.u
        dt = self.lb.dt
        tau = self.lb.tau

        # store original f field
        f = self.lb.f

        n_offset = 1
        for yi in range(n_offset, self.lb.ny-n_offset):
            for xi in range(n_offset, self.lb.nx-n_offset):
                                
                # Skip over cells containing barriers
                if (self.lb.wall[yi, xi]):
                    continue
                else:
                    # Compute the macroscopic density
                    self.lb.rho[yi, xi] = np.sum(self.lb.f[yi,xi,:])

                    # Compute the macroscopic velocities
                    if (self.lb.rho[yi,xi] > 0):
                        self.lb.u[yi,xi,:] = np.dot(self.lb.f[yi,xi,1:], e_[1:]) * (1-(self.lb.rho[yi,xi]-1)+((self.lb.rho[yi,xi]-1)**2))

                    for qi in range(1,9):
                        f_eq = w_[qi] * (1 + 3 * np.dot(u[yi,xi,:], e_[qi]) + 4.5 * (np.dot(u[yi,xi,:], e_[qi]))**2 - 1.5 * np.dot(u[yi,xi,:], u[yi,xi,:]))
                        self.lb.f[yi, xi, qi] =  (1-omega) * f[yi, xi, qi] + omega * f_eq
                        self.lb.f[yi, xi, qi] =  (1-dt/tau) * f[yi, xi, qi] + (dt/tau) * f_eq
                    
                    # Conserve mass
                    self.lb.f[yi, xi, 0]  = self.lb.rho[yi, xi] - np.sum(f[yi, xi, 1:])


    def bounce(self):

        n_offset = 2 
        for yi in range(n_offset, self.lb.ny-n_offset):
            for xi in range(n_offset, self.lb.nx-n_offset):
                
                # If the cell contains a wall
                if (self.lb.wall[yi, xi]):
                    
                    # bounce back
                    self.lb.f[yi+1, xi, 2] = self.lb.f[yi, xi, 4]
                    self.lb.f[yi-1, xi, 4] = self.lb.f[yi, xi, 2]
                    self.lb.f[yi, xi+1, 1] = self.lb.f[yi, xi, 3]
                    self.lb.f[yi, xi-1, 3] = self.lb.f[yi, xi, 1]
                    self.lb.f[yi+1, xi+1, 5] = self.lb.f[yi, xi, 7]
                    self.lb.f[yi+1, xi-1, 6] = self.lb.f[yi, xi, 8]
                    self.lb.f[yi-1, xi+1, 8] = self.lb.f[yi, xi, 6]
                    self.lb.f[yi-1, xi-1, 7] = self.lb.f[yi, xi, 5]
                    
                    # clear f in wall
                    for fi in range(9):
                        self.lb.f[yi, xi, fi] = 0