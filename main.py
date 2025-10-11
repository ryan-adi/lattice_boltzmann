from common_modules import np, plt, time

from src.lbm import *
from src.boundary_condition import * 
from src.postprocessing import *

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
dt = 0.001
tSim = 0.1
nt = int(tSim / dt)
case_name = "channel2D"

# Geometry
grid_size = np.array([64, 16])               # grid dimensions
nx = grid_size[0]                           # grid width
ny = grid_size[1]                           # grid height

# fluid properties
u0 = 0.3
lbm_properties = {
    'viscosity': 0.002,                 # viscosity
    #'omega': 1./(3*viscosity + 0.5),   # relaxation parameter (a function of viscosity)
    'omega': 0.8,                       # source: me
    'c': 1.0,                           # lattice velocity (currently not used)
    'u0': .0 * np.array([1.0, 0.0]),    # initial in-flow velocity (in percentage of c)
    'rho': 1.0,                         # base density
    'dt': dt,                           # time step
    'tau':1.0                           # relaxation time
}

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9


if __name__ == '__main__':

    start = time.time()

    # create png output folder   
    create_output_folder(case_name)

    # initialize LatticeBoltzmann
    lb = LatticeBoltzmann(D, Q)
    lb.init_quantities(lbm_properties)
    lb.init_grid(lx=1.0, ly=1.0, nx=nx, ny=ny)

    # define boundary conditions
    bc_dict = {"top":[], "bottom":[], "left":[], "right":[]}
    bc_dict["left"].append('velocity')
    bc_dict["left"].append(u0*np.array([.5, 0.0]))
    # bc_dict["left"].append('pressure')
    # bc_dict["left"].append(np.array([1.0, 0]))

    bc_dict["right"].append('velocity')
    bc_dict["right"].append(u0*np.array([.5, 0.0]))

    # bc_dict["top"].append('bounce')
    # bc_dict["top"].append(0)

    # bc_dict["bottom"].append('bounce')
    # bc_dict["bottom"].append(0)
    
    ## init objects as a set of wall cells 
    # lb.init_wall([0,ny-1], [nx-1,ny-1])
    # lb.init_wall([0,0], [nx-1,0])
    # lb.init_wall([nx-2,ny-2], [nx-1,ny-1])

    ## simulation loop
    for ti in range(nt+1):

        # get macro quantities
        rho = lb.get_rho()
        ux = lb.get_velocity()[:,:,0]
        uy = lb.get_velocity()[:,:,1]
        ke = lb.get_kinetic_energy()

        # save pngs 
        ds = {
            # "Density":rho, 
              "Velocity X":ux,
              "Velocity Y":uy,
              "Kinetic Energy":ke
              }
    
        # output visualization
        save_pngs(case_name, current_iter=ti, ds=ds)
        save_csvs(case_name, current_iter=ti, ds=ds)

        # LBM update
        lb.update(bc_dict)


    ## generate videos from pngs
    generate_video(case_name, ds=ds, fps=20, remove_pngs=False)

    end = time.time()
    print(f"LBM simulation complete in {end-start} s")