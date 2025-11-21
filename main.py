from common_modules import np, time
import cProfile

from src.lbm import LatticeBoltzmann
from src.initializer import Initializer
from src.obstacle import Obstacle
from src.boundary_condition import * 
from src.postprocess import *
from src.visualization import *

# simulation settings
fps = 10
export_interval = 10
dt = 1.0
tSim = 1000.0
nt = int(tSim / dt)
case_name = "channel2d"

# Geometry
grid_size = np.array([100, 50])            # grid dimensions
nx = grid_size[0]                           # grid width
ny = grid_size[1]                           # grid height

# fluid properties
u0 = 0.2
lbm_properties = {
    'viscosity': 0.002,                 # viscosity
    'c': 0.577,                           # lattice velocity
    'u0': u0 * np.array([1.0, 0.0]),    # initial in-flow velocity (in percentage of c)
    'rho': 1.0,                         # base density
    'dt': dt,                           # time step
    'tau':0.8                           # relaxation time
}

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9

def simulation(export_interval:int) -> None:

    start = time.time()

    # create results folders
    case_dir = os.path.join("output", case_name)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)  
    create_output_folder(case_name)

    # initialize LatticeBoltzmann
    lb = LatticeBoltzmann(D, Q)
    initializer = Initializer(lb)
    initializer.physical_quantities(lbm_properties)
    initializer.grid_quantities(lx=1.0, ly=1.0, nx=nx, ny=ny)
    initializer.micro_velocities(e=[1,2], val=[1.0,1.0])
    initializer.macro_quantities(u0=1.0, rho0=1.0)

    ## create obstacles
    obstacle = Obstacle(lb)
    # obstacle.create_box([nx//3-2,ny//2-10], [nx//3+2,ny//2+10])
    # obstacle.create_circle(radius=5,center=[nx//3,ny//2])

    # define boundary conditions
    bc_dict = {"top":[], "bottom":[], "left":[], "right":[]}
    bc_dict["left"].append('set')
    bc_dict["left"].append(np.array([3.0,1.0]))
    # bc_dict["left"].append('absorb')
    # bc_dict["left"].append([])
    bc_dict["right"].append('absorb')
    bc_dict["right"].append([])
    bc_dict["top"].append('velocity')
    bc_dict["top"].append(np.array([0.0, 0.0]))
    bc_dict["bottom"].append('velocity')
    bc_dict["bottom"].append(np.array([0.0, 0.0]))
    
    # postprocessing
    postprocess = Postprocess(lb, case_name)
    postprocess.create_postprocessing_folder()

    # initialize field by running simulation loop for 10 iterations
    for _ in range(10):
        lb.update(bc_dict)

    ## simulation loop
    for iter in range(nt+1):

        if (iter%export_interval==0):
            print(f"Time = {iter*dt:.2f} s")

        # get macro quantities
        rho = lb.get_rho()
        ux = lb.get_velocity()[:,:,0]
        uy = lb.get_velocity()[:,:,1]
        ke = lb.get_kinetic_energy()
        cu = lb.get_curl()
        scalar = lb.get_scalar()

        # save pngs 
        ds = {
              "Density":rho, 
              "Velocity_X":ux,
              "Velocity_Y":uy,
              "Kinetic_Energy":ke,
              "Curl":cu,
              "Scalar":scalar,
              "F1":lb.f[:,:,1],
              }
    
        # output visualization
        if (iter%export_interval==0): 
            save_pngs(case_name, current_time = iter*dt, 
                      export_iter=iter//export_interval, ds=ds)
            save_csvs(case_name, export_iter=iter//export_interval, ds=ds)
            
            # postprocessing
            postprocess.run(iter)
        
        # LBM update
        lb.update(bc_dict)

        # change f3 
        # if (iter%1==0 & iter<50):
        #     slc1 = slice(ny//2-10, ny//2+10)
        #     slc2 = slice(nx//3-2+iter//1, nx//3+2+iter//1)
        #     lb.f[slc1,slc2,1] += 0.1
        # move obstacles
        # if (iter==10):
        #     obstacle.move([0,1])
        # f = 1.0
        # if (iter%export_interval==0):
        #     direction =  4 * np.floor(f*iter/export_interval) - 2 * np.floor(2*f*iter/export_interval) + 1
        #     displacement = np.array([1, 0]) * direction
        #     obstacle.move(displacement)

    end = time.time()

    ## generate videos from pngs
    generate_video(case_name, ds=ds, fps=fps, remove_pngs=True)

    print(f"LBM simulation complete in {end-start} s")

if __name__ == "__main__":
    print("Starting simulation.")

    ## run code with output
    simulation(export_interval=export_interval)

    ## profile code (visualize with snakeviz lbm.prof)
    #cProfile.run("simulation(export_interval=100)", "lbm.prof")