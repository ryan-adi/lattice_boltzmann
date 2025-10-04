import numpy as np

## constants
c1 = 2.0/3.0
c2 = 1.0/6.0
c3 = 1.0/2.0

class BoundaryCondition():
    def __init__(self, lb):
        self.lb = lb
        
    # ==================== PRESSURE BCs (zou he) ==================== #
    ### left pressure BC
    # def left_pressure_bc(self, u_bc):

    #     first_bc_cell = self.lb.nx
    #     last_bc_cell = (self.lb.ny-2)*self.lb.nx + first_bc_cell
    #     bc_cells = np.arange(first_bc_cell, last_bc_cell, self.lb.nx)

        # self.lb.rho[lx,:] = rho_right[:]
        # self.lb.u[1,lx,:] = u_right[1,:]

        # self.lb.u[0,lx,:] = (g[0,lx,:] + g[3,lx,:] + g[4,lx,:] +
        #              2.0*g[1,lx,:] + 2.0*g[5,lx,:] +
        #              2.0*g[8,lx,:])/self.lb.rho[lx,:] - 1.0

        # g[2,lx,:] = (g[1,lx,:] - cst1*self.lb.rho[lx,:]*self.lb.u[0,lx,:])

        # g[6,lx,:] = (g[5,lx,:] + cst3*(g[3,lx,:] - g[4,lx,:]) -
        #              cst2*self.lb.rho[lx,:]*self.lb.u[0,lx,:] -
        #              cst3*self.lb.rho[lx,:]*self.lb.u[1,lx,:] )

        # g[7,lx,:] = (g[8,lx,:] - cst3*(g[3,lx,:] - g[4,lx,:]) -
        #              cst2*self.lb.rho[lx,:]*self.lb.u[0,lx,:] +
        #              cst3*self.lb.rho[lx,:]*self.lb.u[1,lx,:] )


    # ==================== VELOCITY BCs (zou he) ==================== #
    ### left velocity BC
    def left_velocity_bc(self, u_bc):

        first_bc_cell = self.lb.nx
        last_bc_cell = (self.lb.ny-2)*self.lb.nx + first_bc_cell
        bc_cells = np.arange(first_bc_cell, last_bc_cell, self.lb.nx)

        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells,0] + self.lb.f[bc_cells,2] + self.lb.f[bc_cells,4] 
                        + 2.0*self.lb.f[bc_cells,3] + 2.0*self.lb.f[bc_cells,7] 
                        + 2.0*self.lb.f[bc_cells,6])/(1.0 - self.lb.u[bc_cells,0])
        
        self.lb.f[bc_cells,1] = (self.lb.f[bc_cells,3] + c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0])

        self.lb.f[bc_cells,5] = (self.lb.f[bc_cells,7] 
                        - c3*(self.lb.f[bc_cells,2] - self.lb.f[bc_cells,4]) 
                        + c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] 
                        + c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1])

        self.lb.f[bc_cells,8] = (self.lb.f[bc_cells,6] 
                        + c3*(self.lb.f[bc_cells,2] - self.lb.f[bc_cells,4]) 
                        + c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] 
                        - c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1])
        
        
    ### right velocity BC
    def right_velocity_bc(self, u_bc):

        first_bc_cell = 2*self.lb.nx-1
        last_bc_cell = (self.lb.ny-2)*self.lb.nx + first_bc_cell
        bc_cells = np.arange(first_bc_cell, last_bc_cell, self.lb.nx)
        
        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells,0] + self.lb.f[bc_cells,2] + self.lb.f[bc_cells,4] 
                        + 2.0*self.lb.f[bc_cells,1] + 2.0*self.lb.f[bc_cells,5] 
                        + 2.0*self.lb.f[bc_cells,8])/(1.0 + self.lb.u[bc_cells,0])

        self.lb.f[bc_cells,3] = (self.lb.f[bc_cells,1] - c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0])

        self.lb.f[bc_cells,7] = (self.lb.f[bc_cells,5] 
                        + c3*(self.lb.f[bc_cells,2] - self.lb.f[bc_cells,4]) 
                        - c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] 
                        - c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )

        self.lb.f[bc_cells,6] = (self.lb.f[bc_cells,8] 
                        - c3*(self.lb.f[bc_cells,2] - self.lb.f[bc_cells,4]) 
                        - c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] 
                        + c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )
        

    ### top velocity BC
    def top_velocity_bc(self, u_bc):

        first_bc_cell = (self.lb.ny-1)*self.lb.nx 
        last_bc_cell = first_bc_cell + (self.lb.nx-2)
        bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)

        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells,0] + self.lb.f[bc_cells,1] + self.lb.f[bc_cells,3] +
                    2.0*self.lb.f[bc_cells,2] + 2.0*self.lb.f[bc_cells,5] +
                    2.0*self.lb.f[bc_cells,6])/(1.0 + self.lb.u[bc_cells,1])

        self.lb.f[bc_cells,4] = (self.lb.f[bc_cells,2] - c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1])

        self.lb.f[bc_cells,8] = (self.lb.f[bc_cells,6] - c3*(self.lb.f[bc_cells,1] - self.lb.f[bc_cells,3]) +
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] -
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )

        self.lb.f[bc_cells,7] = (self.lb.f[bc_cells,5] + c3*(self.lb.f[bc_cells,1] - self.lb.f[bc_cells,3]) -
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] -
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )

    ### bottom velocity BC
    def bottom_velocity_bc(self, u_bc):

        first_bc_cell = 1
        last_bc_cell = first_bc_cell + (self.lb.nx-2)
        bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)

        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.rho[bc_cells] = (self.lb.f[bc_cells,0] + self.lb.f[bc_cells,1] + self.lb.f[bc_cells,3] +
                    2.0*self.lb.f[bc_cells,4] + 2.0*self.lb.f[bc_cells,7] +
                    2.0*self.lb.f[bc_cells,8] )/(1.0 - self.lb.u[bc_cells,1])

        self.lb.f[bc_cells,2] = (self.lb.f[bc_cells,4] + c1*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1])

        self.lb.f[bc_cells,5] = (self.lb.f[bc_cells,7] - c3*(self.lb.f[bc_cells,1] - self.lb.f[bc_cells,3]) +
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] +
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )

        self.lb.f[bc_cells,6] = (self.lb.f[bc_cells,8] + c3*(self.lb.f[bc_cells,1] - self.lb.f[bc_cells,3]) -
                    c3*self.lb.rho[bc_cells]*self.lb.u[bc_cells,0] +
                    c2*self.lb.rho[bc_cells]*self.lb.u[bc_cells,1] )
        

    ### ==================== Set Velocity BC ==================== #####
    def set_left_velocity(self, u_bc):
        first_bc_cell = self.lb.nx
        last_bc_cell = (self.lb.ny-2)*self.lb.nx + first_bc_cell
        bc_cells = np.arange(first_bc_cell, last_bc_cell, self.lb.nx)

        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.f[bc_cells,3] = 0.0
        self.lb.f[bc_cells,7] = 0.0
        self.lb.f[bc_cells,6] = 0.0

    def set_right_velocity(self, u_bc):

        first_bc_cell = 2*self.lb.nx-1
        last_bc_cell = (self.lb.ny-2)*self.lb.nx + first_bc_cell
        bc_cells = np.arange(first_bc_cell, last_bc_cell, self.lb.nx)
        
        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.f[bc_cells,3] = 0.0
        self.lb.f[bc_cells,7] = 0.0
        self.lb.f[bc_cells,6] = 0.0

    def set_top_velocity(self, u_bc):

        first_bc_cell = (self.lb.ny-1)*self.lb.nx 
        last_bc_cell = first_bc_cell + (self.lb.nx-2)
        bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)
        
        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.f[bc_cells,4] = 0.0
        self.lb.f[bc_cells,7] = 0.0
        self.lb.f[bc_cells,8] = 0.0

    def set_bottom_velocity(self, u_bc):

        first_bc_cell = 1
        last_bc_cell = first_bc_cell + (self.lb.nx-2)
        bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)
        
        self.lb.u[bc_cells,0] = u_bc[0]
        self.lb.u[bc_cells,1] = u_bc[1]

        self.lb.f[bc_cells,2] = 0.0
        self.lb.f[bc_cells,5] = 0.0
        self.lb.f[bc_cells,6] = 0.0