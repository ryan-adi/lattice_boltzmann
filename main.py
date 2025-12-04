from common_modules import np, time
import cProfile
import argparse

from src.lbm import LatticeBoltzmann
from src.initializer import Initializer
from src.obstacle import Obstacle
from src.xml_reader import XMLReader
from src.boundary_condition import * 
from src.postprocess import *
from src.visualization import *

## LBM PARAMS (CURRENTLY ONLY FOR D=2)
D = 2
Q = 9

def simulation() -> None:
    # simulation timer
    start = time.time()

    # read configuration xml
    xmlreader = XMLReader()
    xmlreader.read("simulation_ctrl.xml")
    ctrl_params = xmlreader.get()

    # simulation configuration
    case_name = ctrl_params["caseName"]
    simsetting = ctrl_params["SimSetting"]
    export_interval = simsetting["export"]
    dt = simsetting["dt"]
    endTime = simsetting["endTime"]
    nt = int(endTime / dt)

    # create results folders
    case_dir = os.path.join("output", case_name)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)  
    create_output_folder(case_name)
    ctrl_name = "simulation_ctrl.xml"
    shutil.copyfile(ctrl_name, os.path.join(case_dir, ctrl_name))
    
    # initialize LatticeBoltzmann
    lb = LatticeBoltzmann(D, Q)
    initializer = Initializer(lb)
    initializer.physical_quantities(ctrl_params["Fluid"])
    initializer.grid_quantities(ctrl_params["Geometry"])
    initializer.field_quantities(ctrl_params["InitialConditions"])

    ## create obstacles
    obstacle = Obstacle(lb)
    if "Obstacle" in ctrl_params.keys():
        obstacle.create_obstacle(ctrl_params["Obstacle"])

    # define boundary conditions
    bc_dict = ctrl_params["BoundaryConditions"]

    # postprocessing
    postprocess = Postprocess(lb, case_name)
    postprocess.create_postprocessing_folder()

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

        # save pngs 
        ds = {
              "Density":rho, 
              "Velocity_X":ux,
              "Velocity_Y":uy,
              "Kinetic_Energy":ke,
              #"Curl":cu,
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
    generate_video(case_name, ds=ds, fps=simsetting["fps"], remove_pngs=simsetting["removePng"])

    print(f"LBM simulation complete in {end-start} s")

if __name__ == "__main__":
    print("Starting simulation.")

    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", action='store_true')
    args = parser.parse_args()

    if args.profile:
        # profile code (visualize with snakeviz lbm.prof)
        cProfile.run("simulation()", "lbm.prof")
    else:
        ## run code with output
        simulation()