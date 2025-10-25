from common_modules import np, plt, time
import cProfile
import snakeviz

from src.lbm import *
from src.boundary_condition import * 
from src.postprocessing import *

# simulation settings
fps = 20
dt = 0.01
tSim = 1
nt = int(tSim / dt)
case_name = "temp2d"

# Geometry
grid_size = np.array([128, 32])              # grid dimensions
nx = grid_size[0]                           # grid width
ny = grid_size[1]                           # grid height

# fluid properties
u0 = 0.3
lbm_properties = {
    'viscosity': 0.002,                 # viscosity
    #'omega': 1./(3*viscosity + 0.5),   # relaxation parameter (a function of viscosity)
    'omega': 0.8,                       # source: me
    'c': 0.577,                           # lattice velocity
    'u0': .0 * np.array([1.0, 0.0]),    # initial in-flow velocity (in percentage of c)
    'rho': 1.0,                         # base density
    'dt': dt,                           # time step
    'tau':1.0                           # relaxation time
}

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9

def simulation(fps:int, postproc:bool) -> None:

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
    bc_dict["left"].append(u0*np.array([1.0, 0.0]))

    # bc_dict["left"].append('pressure')
    # bc_dict["left"].append(np.array([1.0, 0]))

    # bc_dict["right"].append('velocity')
    # bc_dict["right"].append(u0*np.array([.5, 0.0]))

    # bc_dict["top"].append('velocity')
    # bc_dict["top"].append(np.array([0.0, 0.0]))

    # bc_dict["right"].append('velocity')
    # bc_dict["right"].append(np.array([0.0, 0.0]))

    # bc_dict["bottom"].append('velocity')
    # bc_dict["bottom"].append(np.array([0.0, 0.0]))
    
    ## init objects as a set of wall cells 
    # lb.init_wall([nx//3-2,ny//2-2], [nx//3+2,ny//2+2])

    # initialize field by running simulation loop for 10 iterations
    for _ in range(10):
        lb.update(bc_dict)

    ## simulation loop
    for ti in range(nt+1):

        print(f"Time = {ti*dt*1000:.2f} ms")

        # get macro quantities
        rho = lb.get_rho()
        ux = lb.get_velocity()[:,:,0]
        uy = lb.get_velocity()[:,:,1]
        ke = lb.get_kinetic_energy()

        # save pngs 
        ds = {
              "Density":rho, 
              "Velocity X":ux,
              "Velocity Y":uy,
              "Kinetic Energy":ke
              }
    
        # output visualization
        if (postproc):
            save_pngs(case_name, current_iter=ti, ds=ds)
            save_csvs(case_name, current_iter=ti, ds=ds)
        
        # LBM update
        lb.update(bc_dict)

    ## generate videos from pngs
    if (postproc):
        generate_video(case_name, ds=ds, fps=fps, remove_pngs=False)

    end = time.time()
    print(f"LBM simulation complete in {end-start} s")

if __name__ == "__main__":
    print("Starting simulation.")

    ## run code with output 
    simulation(fps=10, postproc=True)

    ## profile code (visualize with snakeviz lbm.prof)
    # cProfile.run("simulation(postproc=False)", "lbm.prof")