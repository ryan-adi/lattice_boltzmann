import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from lbm import *
from boundary_condition import * 

plt.rcParams["figure.figsize"] = (50,3)
plt.rcParams.update({
    'axes.labelsize': 14,      # Font size for axis labels
    'axes.titlesize': 16,      # Font size for axis titles
    'xtick.labelsize': 12,     # Font size for x-axis ticks
    'ytick.labelsize': 12,     # Font size for y-axis ticks
    #'font.family': 'serif',    # Font family (default is 'sans-serif')
    'font.size': 12,           # General font size for text in the plot
})

# simulation settings
fps = 20
nSeconds = 5

# fluid properties
viscosity = 0.002                       # viscosity
omega = 1./(3*viscosity + 0.5)          # relaxation parameter (a function of viscosity)
omega = 0.8                             # source: me
u0 = 0.0 * np.array([1.0, 0.0])         # initial in-flow velocity
c = 1.0                                 # lattice velocity (currently not used)
#rho_base = 1000                        # base density

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9
# for D2Q9
e_ = np.array([np.array([0.0, 0.0]),
                np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
for e_i in e_[1:]: # normalize vectors
    e_i /= np.sqrt(np.dot(e_i, e_i))
w_ = np.array([4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.])

# Geometry
grid_size = np.array([64, 16])              # grid dimensions
width = grid_size[0]                        # grid width
height = grid_size[1]                       # grid height

# Equilbrium distributions
f = np.zeros((height*width, Q))

# Walls
wall = np.zeros(height*width)

# Macroscopic density and velocity
rho = np.zeros(height*width)    # density
u = np.zeros((height*width, D))   # velocity in 3D
speed2 = np.zeros(height*width) # squared velocity

def initialize(xtop, ytop, yheight, u0=u0):
    xcoord = 0
    ycoord = 0
    
    count = 0
    for i in range(height*width):
        
        # init equilibrium distribution
        for qi in range(Q):
            f[i, qi] = w_[qi] * (1 + 3 * np.dot(u0, e_[qi]) + 4.5 * (np.dot(u0, e_[qi]))**2 - 1.5 * np.dot(u0, u0)) 
        
        # macroscopic values
        rho[i] = np.sum(f[i,:])
        u[i,:] = np.dot(f[i,1:], e_[1:]) * (1-(rho[i]-1)+((rho[i]-1)**2))
        
        # # wall init
        # if (xcoord==xtop):
        #     if (ycoord >= ytop):
        #         if (ycoord < (ytop+yheight)):
        #             count += 1
        #             wall[ycoord*width + xcoord] = 1
        
        # xcoord = (xcoord+1) if xcoord<(width-1) else 0
        # ycoord = ycoord if (xcoord != 0) else (ycoord + 1)


def print_csv(f): 
    pass

def set_colorbar(im, label):
    cbar = plt.colorbar(im)
    plt.gca().invert_yaxis()
    cbar.set_label(label, rotation=270)


if __name__ == '__main__':
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(20,5))

    # Initialize the Walls (occurs in previous section)
    initialize(xtop=25, ytop=11, yheight=10)

    # start first 10 iteration (for testing purposes)
    for i in range(10):
        ## Boundary Conditions
        # left_velocity_bc(width, height, f, u, speed2, rho, np.array([0.0, 0.0]))

        ## LBM operations
        stream(width, height, f)
        bounce(width, height, f, wall)
        collide(width, height, w_, e_, f, rho, u, speed2, wall, omega)
        # left_wall_velocity_bc(width, height, u, u0, rho, f)

    a = speed2 #u[:,0]
    # a = np.zeros(a.shape)
    # a[(width-1)*height] = 100
    # a[16] = 100
    im = plt.imshow(a.reshape(height,width))
    set_colorbar(im, label="speed2")

    # Animation function (stream, bounce, collide, and update heatmap)
    def animate_func(i):
        
        ## LBM operations
        stream(width, height, f)
        bounce(width, height, f, wall)
        collide(width, height, w_, e_, f, rho, u, speed2, wall, omega)

        ## Boundary Conditions
        left_velocity_bc(width, height, f, u, speed2, rho, np.array([0.1, 0.0]))
        right_velocity_bc(width, height, f, u, speed2, rho, np.array([0.1, 0.0]))
        # top_velocity_bc(width, height, f, u, speed2, rho, np.array([0.0, 0.0]) )
        # bottom_velocity_bc(width, height, f, u, speed2, rho, np.array([0.0, 0.0]) )

        im.set_data(a.reshape(height, width))
        plt.xlabel("x coordinate")
        plt.ylabel("y coordinate")
        return [im]

    # Animation object
    anim = animation.FuncAnimation(
                                fig, 
                                animate_func, 
                                frames = nSeconds * fps,
                                interval = 1000 / fps, # in ms
                                )

    print('Done with calculations!')

    # Generate an mp4 video of the animation
    print("Generate video")
    vid_type = "gif"
    if (vid_type=="mp4"):
        file_name = r"./animation4.mp4" 
    # writervideo = animation.FFMpegWriter(fps=fps) 
    else:
        file_name = r"./animation4.gif" 
        writervideo = animation.PillowWriter(fps=600)
    anim.save(file_name, writer=writervideo)
    print("Done generating video")