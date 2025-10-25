from common_modules import np, plt, time
import cProfile

from src.lbm import LatticeBoltzmann
from src.initializer import Initializer
from src.obstacle import Obstacle
from src.boundary_condition import * 
from src.postprocessing import *

# simulation settings
fps = 20
export = 50
dt = 1.0
tSim = 5000.0
nt = int(tSim / dt)
case_name = "temp2d"

# Geometry
grid_size = np.array([750, 250])            # grid dimensions
nx = grid_size[0]                           # grid width
ny = grid_size[1]                           # grid height

# fluid properties
u0 = 0.3
lbm_properties = {
    'viscosity': 0.002,                 # viscosity
    'c': 0.577,                           # lattice velocity
    'u0': .0 * np.array([1.0, 0.0]),    # initial in-flow velocity (in percentage of c)
    'rho': 1.0,                         # base density
    'dt': dt,                           # time step
    'tau':0.53                           # relaxation time
}

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9

def simulation(export_interval:int) -> None:

    start = time.time()

    # create png output folder   
    create_output_folder(case_name)

    # initialize LatticeBoltzmann
    lb = LatticeBoltzmann(D, Q)
    initializer = Initializer(lb)
    initializer.physical_quantities(lbm_properties)
    initializer.grid_quantities(lx=1.0, ly=1.0, nx=nx, ny=ny)
    initializer.micro_velocities(e=[1,2], val=[3.0,1.0])
    initializer.macro_quantities()

    ## create obstacles
    obstacle = Obstacle(lb)
    obstacle.create_box([2*nx//3-2,ny//2-2], [2*nx//3+2,ny//2+2])
    obstacle.create_circle(radius=5,center=[nx//3,ny//2])

    # define boundary conditions
    bc_dict = {"top":[], "bottom":[], "left":[], "right":[]}
    # bc_dict["left"].append('velocity')
    # bc_dict["left"].append(u0*np.array([1.0, 0.0]))

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
    
   
    # initialize field by running simulation loop for 10 iterations
    for _ in range(10):
        lb.update(bc_dict)

    ## simulation loop
    for ti in range(nt+1):

        print(f"Time = {ti*dt:.2f} s")

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
        if (ti%export_interval==0): 
            save_pngs(case_name, current_time = ti*dt, 
                      export_iter=ti//export_interval, ds=ds)
            save_csvs(case_name, current_time=ti*dt, 
                      export_iter=ti//export_interval, ds=ds)
        
        # LBM update
        lb.update(bc_dict)

    end = time.time()

    ## generate videos from pngs
    generate_video(case_name, ds=ds, fps=fps, remove_pngs=True)

    print(f"LBM simulation complete in {end-start} s")

if __name__ == "__main__":
    print("Starting simulation.")

    ## run code with output
    simulation(export_interval=export)

    ## profile code (visualize with snakeviz lbm.prof)
    #cProfile.run("simulation(export_interval=10)", "lbm.prof")