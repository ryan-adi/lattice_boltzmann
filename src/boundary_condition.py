
## constants
c1 = 2.0/3.0
c2 = 1.0/6.0
c3 = 1.0/2.0

class BoundaryCondition():
    def __init__(self, lb):
        self.lb = lb
        
    # ==================== PRESSURE BCs (zou he) ==================== #
    ### left pressure BC
    def pressure_left_bc(self, val):
        """
        Pressure Boundary condition for left boundary

        :param list val: [rho, u[1]]
        :return: None
        """

        bc_cells = (slice(1,self.lb.ny-1),0)

        self.lb.rho[bc_cells] = val[0]
        self.lb.u[bc_cells][:,1] = val[1]

        self.lb.u[bc_cells][:,0] = (self.lb.f[bc_cells][:,0] + self.lb.f[bc_cells][:,3] + self.lb.f[bc_cells][:,4] +
                     2.0*self.lb.f[bc_cells][:,1] + 2.0*self.lb.f[bc_cells][:,5] +
                     2.0*self.lb.f[bc_cells][:,8])/self.lb.rho[bc_cells] - 1.0

        self.lb.f[bc_cells][:,2] = (self.lb.f[bc_cells][:,1] - c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0])

        self.lb.f[bc_cells][:,6] = (self.lb.f[bc_cells][:,5] + c3*(self.lb.f[bc_cells][:,3] - self.lb.f[bc_cells][:,4]) -
                     c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] -
                     c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )

        self.lb.f[bc_cells][:,7] = (self.lb.f[bc_cells][:,8] - c3*(self.lb.f[bc_cells][:,3] - self.lb.f[bc_cells][:,4]) -
                     c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] +
                     c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1])


    # ==================== VELOCITY BCs (zou he) ==================== #
    ### left velocity BC
    def velocity_left_bc(self, u_bc):
        """
        Velocity boundary condition for left boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(1,self.lb.ny-1),0)

        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells][:,0] + self.lb.f[bc_cells][:,2] + self.lb.f[bc_cells][:,4] 
                        + 2.0*self.lb.f[bc_cells][:,3] + 2.0*self.lb.f[bc_cells][:,7] 
                        + 2.0*self.lb.f[bc_cells][:,6])/(1.0 - self.lb.u[bc_cells][:,0])
        
        self.lb.f[bc_cells][:,1] = (self.lb.f[bc_cells][:,3] + c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0])

        self.lb.f[bc_cells][:,5] = (self.lb.f[bc_cells][:,7] 
                        - c3*(self.lb.f[bc_cells][:,2] - self.lb.f[bc_cells][:,4]) 
                        + c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] 
                        + c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1])

        self.lb.f[bc_cells][:,8] = (self.lb.f[bc_cells][:,6] 
                        + c3*(self.lb.f[bc_cells][:,2] - self.lb.f[bc_cells][:,4]) 
                        + c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] 
                        - c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1])
        
        
    ### right velocity BC
    def velocity_right_bc(self, u_bc):
        """
        Velocity boundary condition for right boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(1,self.lb.ny-1), 1) #self.lb.nx-1)
        
        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells][:,0] + self.lb.f[bc_cells][:,2] + self.lb.f[bc_cells][:,4] 
                        + 2.0*self.lb.f[bc_cells][:,1] + 2.0*self.lb.f[bc_cells][:,5] 
                        + 2.0*self.lb.f[bc_cells][:,8])/(1.0 + self.lb.u[bc_cells][:,0])

        self.lb.f[bc_cells][:,3] = (self.lb.f[bc_cells][:,1] - c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0])

        self.lb.f[bc_cells][:,7] = (self.lb.f[bc_cells][:,5] 
                        + c3*(self.lb.f[bc_cells][:,2] - self.lb.f[bc_cells][:,4]) 
                        - c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] 
                        - c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )

        self.lb.f[bc_cells][:,6] = (self.lb.f[bc_cells][:,8] 
                        - c3*(self.lb.f[bc_cells][:,2] - self.lb.f[bc_cells][:,4]) 
                        - c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] 
                        + c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )
        

    ### top velocity BC
    def velocity_top_bc(self, u_bc):
        """
        Velocity boundary condition for top boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (self.lb.ny-1, slice(1, self.lb.nx-1) )

        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells][:,0] + self.lb.f[bc_cells][:,1] + self.lb.f[bc_cells][:,3] +
                    2.0*self.lb.f[bc_cells][:,2] + 2.0*self.lb.f[bc_cells][:,5] +
                    2.0*self.lb.f[bc_cells][:,6])/(1.0 + self.lb.u[bc_cells][:,1])

        self.lb.f[bc_cells][:,4] = (self.lb.f[bc_cells][:,2] - c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1])

        self.lb.f[bc_cells][:,8] = (self.lb.f[bc_cells][:,6] - c3*(self.lb.f[bc_cells][:,1] - self.lb.f[bc_cells][:,3]) +
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] -
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )

        self.lb.f[bc_cells][:,7] = (self.lb.f[bc_cells][:,5] + c3*(self.lb.f[bc_cells][:,1] - self.lb.f[bc_cells][:,3]) -
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] -
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )

    ### bottom velocity BC
    def velocity_bottom_bc(self, u_bc):
        """
        Velocity boundary condition for bottom boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (1, slice(1, self.lb.nx-1))

        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells][:,0] + self.lb.f[bc_cells][:,1] + self.lb.f[bc_cells][:,3] +
                    2.0*self.lb.f[bc_cells][:,4] + 2.0*self.lb.f[bc_cells][:,7] +
                    2.0*self.lb.f[bc_cells][:,8] )/(1.0 - self.lb.u[bc_cells][:,1])

        self.lb.f[bc_cells][:,2] = (self.lb.f[bc_cells][:,4] + c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1])

        self.lb.f[bc_cells][:,5] = (self.lb.f[bc_cells][:,7] - c3*(self.lb.f[bc_cells][:,1] - self.lb.f[bc_cells][:,3]) +
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] +
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )

        self.lb.f[bc_cells][:,6] = (self.lb.f[bc_cells][:,8] + c3*(self.lb.f[bc_cells][:,1] - self.lb.f[bc_cells][:,3]) -
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,0] +
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells][:,1] )
        
    ### ==================== Bounce-Back BC ==================== #####    

    def bounce_back_left_bc(self):

        xi = 0

        for yi in range(0,self.lb.ny):
            self.lb.f[yi, xi+1, 1] = self.lb.f[yi, xi, 3]
            if (yi<self.lb.ny-1):
                self.lb.f[yi+1, xi+1, 5] = self.lb.f[yi, xi, 7]
            if (yi>0):
                self.lb.f[yi-1, xi+1, 8] = self.lb.f[yi, xi, 6]

            # clear f in wall
            for fi in range(9):
                self.lb.f[yi, xi, fi] = 0  

    def bounce_back_right_bc(self):

        xi = self.lb.nx-1

        for yi in range(0,self.lb.ny):
            self.lb.f[yi, xi-1, 3] = self.lb.f[yi, xi, 1]
            if (yi>0):
                self.lb.f[yi-1, xi-1, 7] = self.lb.f[yi, xi, 5]
            if (yi<self.lb.ny-1):
                self.lb.f[yi+1, xi-1, 6] = self.lb.f[yi, xi, 8]

            # clear f in wall
            for fi in range(9):
                self.lb.f[yi, xi, fi] = 0  

    def bounce_back_top_bc(self):

        yi = self.lb.ny-1

        for xi in range(0, self.lb.nx):
            self.lb.f[yi-1, xi, 4] = self.lb.f[yi, xi, 2]
            if (xi>0) :
                self.lb.f[yi-1, xi-1, 7] = self.lb.f[yi, xi, 5]
            if (xi<self.lb.nx-1):
                self.lb.f[yi-1, xi+1, 8] = self.lb.f[yi, xi, 6]

            # clear f in wall
            for fi in range(9):
                self.lb.f[yi, xi, fi] = 0   

    def bounce_back_bottom_bc(self):

        yi = 0

        for xi in range(0, self.lb.nx):
            self.lb.f[yi+1, xi, 2] = self.lb.f[yi, xi, 4]
            if (xi<self.lb.nx-1) :
                self.lb.f[yi+1, xi+1, 5] = self.lb.f[yi, xi, 7]
            if (xi>0) :
                self.lb.f[yi+1, xi-1, 6] = self.lb.f[yi, xi, 8]

            # clear f in wall
            for fi in range(9):
                self.lb.f[yi, xi, fi] = 0  
      
        
    ### ==================== Set Velocity BC ==================== #####
    def set_velocity_left_bc(self, u_bc):
        """
        Set Velocity and f at left boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(1,self.lb.ny-1),0)

        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.f[bc_cells][:,3] = 0.0
        self.lb.f[bc_cells][:,7] = 0.0
        self.lb.f[bc_cells][:,6] = 0.0

    def set_velocity_right_bc(self, u_bc):
        """
        Set Velocity and f at right boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(1,self.lb.ny-1), self.lb.nx-1)
        
        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.f[bc_cells][:,3] = 0.0
        self.lb.f[bc_cells][:,7] = 0.0
        self.lb.f[bc_cells][:,6] = 0.0

    def set_velocity_top_bc(self, u_bc):
        """
        Set Velocity and f at top boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (self.lb.ny-1, slice(1, self.lb.nx-1) )
        
        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.f[bc_cells][:,4] = 0.0
        self.lb.f[bc_cells][:,7] = 0.0
        self.lb.f[bc_cells][:,8] = 0.0

    def set_velocity_bottom_bc(self, u_bc):
        """
        Set Velocity and f at bottom boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (1, slice(1, self.lb.nx-1))
        
        self.lb.u[bc_cells][:,0] = u_bc[0]
        self.lb.u[bc_cells][:,1] = u_bc[1]

        self.lb.f[bc_cells][:,2] = 0.0
        self.lb.f[bc_cells][:,5] = 0.0
        self.lb.f[bc_cells][:,6] = 0.0