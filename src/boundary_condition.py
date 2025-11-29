from common_modules import np

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

        ny = self.lb.ny
        ny0 = 0
        ny1 = ny

        bc_cells = (slice(ny0,ny1),0)

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

        ny = self.lb.ny
        ny0 = 0
        ny1 = ny

        bc_cells = (slice(ny0, ny1),0)

        # self.lb.u[bc_cells][:,0] = u_bc[0]
        # if poiseuille
        # for j in range(ny0, ny1):
        #     v_poi = u_bc[0] * (1- ((j-ny/2)/(ny-2))**2)
        #     self.lb.u[j,0,0] = v_poi

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

        bc_cells = (slice(1,self.lb.ny-1), self.lb.nx-1)
        
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

        bc_cells = (self.lb.ny-1, slice(1, self.lb.nx-1))

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

        bc_cells = (0, slice(1, self.lb.nx-1))

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
        
    ### ==================== Set Microscopic Velocity BC ==================== #####
    def set_left_bc(self, f_bc):
        """
        Set Velocity and f at left boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(0,self.lb.ny),0)

        idx1 = 1
        if (f_bc[0]<0.0):
            idx1 = 3
        idx2 = 2
        if (f_bc[1]<0.0):
            idx2 = 4

        self.lb.f[bc_cells][:,idx1] = abs(f_bc[0])
        self.lb.f[bc_cells][:,idx2] = abs(f_bc[1])

    def set_right_bc(self, f_bc):
        """
        Set Velocity and f at right boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (slice(0,self.lb.ny), self.lb.nx-1)
        
        idx1 = 1
        if (f_bc[0]<0.0):
            idx1 = 3
        idx2 = 2
        if (f_bc[1]<0.0):
            idx2 = 4

        self.lb.f[bc_cells][:,idx1] = abs(f_bc[0])
        self.lb.f[bc_cells][:,idx2] = abs(f_bc[1])

    def set_top_bc(self, f_bc):
        """
        Set Velocity and f at top boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (self.lb.ny-1, slice(0, self.lb.nx) )
        
        idx1 = 1
        if (f_bc[0]<0.0):
            idx1 = 3
        idx2 = 2
        if (f_bc[1]<0.0):
            idx2 = 4

        self.lb.f[bc_cells][:,idx1] = abs(f_bc[0])
        self.lb.f[bc_cells][:,idx2] = abs(f_bc[1])

    def set_bottom_bc(self, f_bc):
        """
        Set Velocity and f at bottom boundary

        :param list val: [u[0], u[1]]
        :return: None
        """

        bc_cells = (1, slice(0, self.lb.nx))
        
        idx1 = 1
        if (f_bc[0]<0.0):
            idx1 = 3
        idx2 = 2
        if (f_bc[1]<0.0):
            idx2 = 4

        self.lb.f[bc_cells][:,idx1] = abs(f_bc[0])
        self.lb.f[bc_cells][:,idx2] = abs(f_bc[1])


    ### ==================== Absorbing BC ==================== #####
    def absorb_left_bc(self):
        self.lb.f[:,0,[1,5,8]] = self.lb.f[:,1,[1,5,8]]  
        
    def absorb_right_bc(self):
        self.lb.f[:,-1,[3,6,7]] = self.lb.f[:,-2,[3,6,7]] 

    def absorb_top_bc(self):
        self.lb.f[-1,:,[4,7,8]] = self.lb.f[-2,:,[4,7,8]] 

    def absorb_bottom_bc(self):
        self.lb.f[0,:,[2,5,6]] = self.lb.f[1,:,[2,5,6]]


    ## ==================== Corner BC ==================== #####
    def corner_bc(self):
        # bottom corners
        for idx1, idx2 in [[0,0],[0,-1]]:
            self.lb.u[idx1,idx2,:] = 0.0
            self.lb.rho[idx1,idx2] = (self.lb.f[idx1,idx2,0] + self.lb.f[idx1,idx2,1] + self.lb.f[idx1,idx2,3] +
                    2.0*self.lb.f[idx1,idx2,4] + 2.0*self.lb.f[idx1,idx2,7] +
                    2.0*self.lb.f[idx1,idx2,8] )/(1.0 - self.lb.u[idx1,idx2,1])
            
            self.lb.f[idx1, idx2,2] = (self.lb.f[idx1, idx2,4] + c1*self.lb.rho[idx1, idx2]*self.lb.u[idx1, idx2,1])

            self.lb.f[idx1, idx2,5] = (self.lb.f[idx1, idx2,7] - c3*(self.lb.f[idx1, idx2,1] - self.lb.f[idx1, idx2,3]) +
                        c3*self.lb.rho[idx1, idx2]*self.lb.u[idx1, idx2,0] +
                        c2*self.lb.rho[idx1, idx2]*self.lb.u[idx1, idx2,1] )

            self.lb.f[idx1, idx2,6] = (self.lb.f[idx1, idx2,8] + c3*(self.lb.f[idx1, idx2,1] - self.lb.f[idx1, idx2,3]) -
                        c3*self.lb.rho[idx1, idx2]*self.lb.u[idx1, idx2,0] +
                        c2*self.lb.rho[idx1, idx2]*self.lb.u[idx1, idx2,1] )
            
        # top corners    
        for idx1, idx2 in [[-1,0],[-1,-1]]:
            self.lb.u[idx1,idx2,:] = 0.0
            self.lb.rho[idx1,idx2] = (self.lb.f[idx1,idx2,0] + self.lb.f[idx1,idx2,1] + self.lb.f[idx1,idx2,3] +
                    2.0*self.lb.f[idx1,idx2,2] + 2.0*self.lb.f[idx1,idx2,5] +
                    2.0*self.lb.f[idx1,idx2,6])/(1.0 + self.lb.u[idx1,idx2,1])
            
            self.lb.f[idx1,idx2,4] = (self.lb.f[idx1,idx2,2] - c1*self.lb.rho[idx1,idx2]*self.lb.u[idx1,idx2,1])

            self.lb.f[idx1,idx2,8] = (self.lb.f[idx1,idx2,6] - c3*(self.lb.f[idx1,idx2,1] - self.lb.f[idx1,idx2,3]) +
                        c3*self.lb.rho[idx1,idx2]*self.lb.u[idx1,idx2,0] -
                        c2*self.lb.rho[idx1,idx2]*self.lb.u[idx1,idx2,1] )

            self.lb.f[idx1,idx2,7] = (self.lb.f[idx1,idx2,5] + c3*(self.lb.f[idx1,idx2,1] - self.lb.f[idx1,idx2,3]) -
                        c3*self.lb.rho[idx1,idx2]*self.lb.u[idx1,idx2,0] -
                        c2*self.lb.rho[idx1,idx2]*self.lb.u[idx1,idx2,1] )