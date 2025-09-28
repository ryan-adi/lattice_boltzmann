import numpy as np

## constants
c1 = 2.0/3.0
c2 = 1.0/6.0
c3 = 1.0/2.0

# ---------------- PRESSURE BCs ---------------- #
### left pressure BC
# def left_pressure_bc(width, height, f, u, rho, u_bc):

#     first_bc_cell = width
#     last_bc_cell = (height-2)*width + first_bc_cell
#     bc_cells = np.arange(first_bc_cell, last_bc_cell, width)

    # rho[lx,:] = rho_right[:]
    # u[1,lx,:] = u_right[1,:]

    # u[0,lx,:] = (g[0,lx,:] + g[3,lx,:] + g[4,lx,:] +
    #              2.0*g[1,lx,:] + 2.0*g[5,lx,:] +
    #              2.0*g[8,lx,:])/rho[lx,:] - 1.0

    # g[2,lx,:] = (g[1,lx,:] - cst1*rho[lx,:]*u[0,lx,:])

    # g[6,lx,:] = (g[5,lx,:] + cst3*(g[3,lx,:] - g[4,lx,:]) -
    #              cst2*rho[lx,:]*u[0,lx,:] -
    #              cst3*rho[lx,:]*u[1,lx,:] )

    # g[7,lx,:] = (g[8,lx,:] - cst3*(g[3,lx,:] - g[4,lx,:]) -
    #              cst2*rho[lx,:]*u[0,lx,:] +
    #              cst3*rho[lx,:]*u[1,lx,:] )


# ---------------- VELOCITY BCs ---------------- #
### left velocity BC
def left_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = width
    last_bc_cell = (height-2)*width + first_bc_cell
    bc_cells = np.arange(first_bc_cell, last_bc_cell, width)

    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    rho[bc_cells] = (f[bc_cells,0] + f[bc_cells,2] + f[bc_cells,4] 
                     + 2.0*f[bc_cells,3] + 2.0*f[bc_cells,7] 
                     + 2.0*f[bc_cells,6])/(1.0 - u[bc_cells,0])
    
    f[bc_cells,1] = (f[bc_cells,3] + c1*rho[bc_cells]*u[bc_cells,0])

    f[bc_cells,5] = (f[bc_cells,7] 
                     - c3*(f[bc_cells,2] - f[bc_cells,4]) 
                     + c2*rho[bc_cells]*u[bc_cells,0] 
                     + c3*rho[bc_cells]*u[bc_cells,1])

    f[bc_cells,8] = (f[bc_cells,6] 
                     + c3*(f[bc_cells,2] - f[bc_cells,4]) 
                     + c2*rho[bc_cells]*u[bc_cells,0] 
                     - c3*rho[bc_cells]*u[bc_cells,1])
    
    
### right velocity BC
def right_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = 2*width-1
    last_bc_cell = (height-2)*width + first_bc_cell
    bc_cells = np.arange(first_bc_cell, last_bc_cell, width)
    
    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    rho[bc_cells] = (f[bc_cells,0] + f[bc_cells,2] + f[bc_cells,4] 
                     + 2.0*f[bc_cells,1] + 2.0*f[bc_cells,5] 
                     + 2.0*f[bc_cells,8])/(1.0 + u[bc_cells,0])

    f[bc_cells,3] = (f[bc_cells,1] - c1*rho[bc_cells]*u[bc_cells,0])

    f[bc_cells,7] = (f[bc_cells,5] 
                     + c3*(f[bc_cells,2] - f[bc_cells,4]) 
                     - c2*rho[bc_cells]*u[bc_cells,0] 
                     - c3*rho[bc_cells]*u[bc_cells,1] )

    f[bc_cells,6] = (f[bc_cells,8] 
                     - c3*(f[bc_cells,2] - f[bc_cells,4]) 
                     - c2*rho[bc_cells]*u[bc_cells,0] 
                     + c3*rho[bc_cells]*u[bc_cells,1] )
    

### top velocity BC
def top_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = (height-1)*width 
    last_bc_cell = first_bc_cell + (width-2)
    bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)

    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    rho[bc_cells] = (f[bc_cells,0] + f[bc_cells,1] + f[bc_cells,3] +
                2.0*f[bc_cells,2] + 2.0*f[bc_cells,5] +
                2.0*f[bc_cells,6])/(1.0 + u[bc_cells,1])

    f[bc_cells,4] = (f[bc_cells,2] - c1*rho[bc_cells]*u[bc_cells,1])

    f[bc_cells,8] = (f[bc_cells,6] - c3*(f[bc_cells,1] - f[bc_cells,3]) +
                 c3*rho[bc_cells]*u[bc_cells,0] -
                 c2*rho[bc_cells]*u[bc_cells,1] )

    f[bc_cells,7] = (f[bc_cells,5] + c3*(f[bc_cells,1] - f[bc_cells,3]) -
                 c3*rho[bc_cells]*u[bc_cells,0] -
                 c2*rho[bc_cells]*u[bc_cells,1] )

### bottom velocity BC
def bottom_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = 1
    last_bc_cell = first_bc_cell + (width-2)
    bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)

    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    rho[bc_cells] = (f[bc_cells,0] + f[bc_cells,1] + f[bc_cells,3] +
                2.0*f[bc_cells,4] + 2.0*f[bc_cells,7] +
                2.0*f[bc_cells,8] )/(1.0 - u[bc_cells,1])

    f[bc_cells,2] = (f[bc_cells,4] + c1*rho[bc_cells]*u[bc_cells,1])

    f[bc_cells,5] = (f[bc_cells,7] - c3*(f[bc_cells,1] - f[bc_cells,3]) +
                c3*rho[bc_cells]*u[bc_cells,0] +
                c2*rho[bc_cells]*u[bc_cells,1] )

    f[bc_cells,6] = (f[bc_cells,8] + c3*(f[bc_cells,1] - f[bc_cells,3]) -
                c3*rho[bc_cells]*u[bc_cells,0] +
                c2*rho[bc_cells]*u[bc_cells,1] )
    

### Set Velocity BC #####
### set right velocity BC
def set_right_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = 2*width-1
    last_bc_cell = (height-2)*width + first_bc_cell
    bc_cells = np.arange(first_bc_cell, last_bc_cell, width)
    
    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    f[bc_cells,3] = 0.0
    f[bc_cells,7] = 0.0
    f[bc_cells,6] = 0.0

def set_top_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = (height-1)*width 
    last_bc_cell = first_bc_cell + (width-2)
    bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)
    
    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    f[bc_cells,4] = 0.0
    f[bc_cells,7] = 0.0
    f[bc_cells,8] = 0.0

def set_bottom_velocity_bc(width, height, f, u, rho, u_bc):

    first_bc_cell = 1
    last_bc_cell = first_bc_cell + (width-2)
    bc_cells = np.arange(first_bc_cell, last_bc_cell, 1)
    
    u[bc_cells,0] = u_bc[0]
    u[bc_cells,1] = u_bc[1]

    f[bc_cells,2] = 0.0
    f[bc_cells,5] = 0.0
    f[bc_cells,6] = 0.0