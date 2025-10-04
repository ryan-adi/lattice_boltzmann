import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc

from src.lbm import *
from src.boundary_condition import * 

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

# Geometry
grid_size = np.array([32, 5])               # grid dimensions
nx = grid_size[0]                           # grid width
ny = grid_size[1]                           # grid height

# fluid properties
physical_properties = {
    'viscosity': 0.002,                 # viscosity
    #'omega': 1./(3*viscosity + 0.5),   # relaxation parameter (a function of viscosity)
    'omega': 0.8,                       # source: me
    'c': 1.0,                           # lattice velocity (currently not used)
    'u0': .5 * np.array([1.0, 0.0]),    # initial in-flow velocity
    #rho_base = 1000                    # base density
}

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9

def set_colorbar(im, label):
    cbar = plt.colorbar(im)
    plt.gca().invert_yaxis()
    cbar.set_label(label, rotation=270)

if __name__ == '__main__':
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(20,5))

    # initialize LatticeBoltzmann
    lb = LatticeBoltzmann(D, Q)
    lb.init_physical_properties(physical_properties)
    lb.init_grid(1.0, 1.0, nx, ny)
    #lb.init_wall([nx-2,ny-2], [nx-1,ny-1])

    ##start first nSeconds iteration (for testing purposes)
    for i in range(nSeconds):

        # time step
        lb.update()
        
        # get macro quantities
        rho = lb.get_rho()
        u = lb.get_velocity()
        ke = lb.get_kinetic_energy()

        # boundary conditions
        # left_velocity_bc(nx, ny, f, u, rho, np.array([0.1, 0.0]))
        # right_velocity_bc(nx, ny, f, u, rho, np.array([0.0, 0.0]))
        # kinetic_energy_update(u, kinetic_energy)

    a = lb.get_kinetic_energy()
    # a = np.zeros(a.shape)
    # a[(nx-1)*ny] = 100
    # a[16] = 100
    im = plt.imshow(a.reshape(ny,nx), vmin=0.0, vmax=0.05, cmap='viridis')
    set_colorbar(im, label="Kinetic Energy")

    # Animation function (stream, bounce, collide, and update heatmap)
    def animate_func(i):
        
        ## LBM update
        #lb.stream()
        #lb.bounce()
        #lb.collide()

        ## Boundary Conditions
        # left_velocity_bc(nx, ny, f, u, rho, np.array([0.1, 0.0]))
        # right_velocity_bc(nx, ny, f, u, rho, np.array([0.1, 0.0]))
        # set_top_velocity_bc(nx, ny, f, u, rho, np.array([0.0, 0.0]))
        # set_bottom_velocity_bc(nx, ny, f, u, rho, np.array([0.0, 0.0]))
        
        ## Derived variable update
        a = lb.get_kinetic_energy()

        ## plot settings
        im.set_data(a.reshape(ny, nx))
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