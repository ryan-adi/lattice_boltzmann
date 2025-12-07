from common_modules import *
import cProfile
import argparse

from src.lbm import LatticeBoltzmann
from src.particles import Particles
from src.obstacle import Obstacle
from src.multiphase import Multiphase
from src.xml_reader import XMLReader
from src.boundary_condition import BoundaryCondition
from src.postprocess import Postprocess
from src.visualization import Visualization

## LBM PARAMS
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
    export_interval = ctrl_params["SimSetting"]["export"]
    dt = ctrl_params["SimSetting"]["dt"]
    endTime = ctrl_params["SimSetting"]["endTime"]
    nt = int(endTime / dt)

    # create results folders
    case_dir = os.path.join("output", case_name)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)  
    ctrl_name = "simulation_ctrl.xml"
    shutil.copyfile(ctrl_name, os.path.join(case_dir, ctrl_name))
    
    # initialize LatticeBoltzmann
    lbm = Multiphase(D, Q)
    lbm.init_physical_quantities(ctrl_params["Fluid"])
    lbm.init_multiphase_components(ctrl_params["Multiphase"])
    lbm.init_grid_quantities(ctrl_params["Geometry"])
    lbm.init_field_quantities(ctrl_params["InitialConditions"])

    # define boundary conditions
    bc_dict = ctrl_params["BoundaryConditions"]

    ## create obstacles
    obstacle = Obstacle(lbm)
    if "Obstacle" in ctrl_params.keys():
        obstacle.create_obstacle(ctrl_params["Obstacle"])
 
    ## Particles
    particles = Particles(x_bound=[0,lbm.nx], y_bound=[0,lbm.ny])
    if "Particles" in ctrl_params.keys():
        particles.initialize(ctrl_params["Particles"])

    # Visualization
    ctrl_params["Visualization"]["CaseName"] = case_name
    visualization = Visualization(lbm, ctrl_params["Visualization"])
    visualization.create_output_folder()

    # postprocessing
    postprocess = Postprocess(lbm, case_name)
    postprocess.create_postprocessing_folder()

    ## simulation loop
    for iter in range(nt+1):
        if (iter%export_interval==0):
            print(f"Time = {iter*dt:.2f} s")

        # save pngs 
        ds = {
              "Density":lbm.get_rho(), 
              "Velocity_X":lbm.get_velocity()[:,:,0],
              "Velocity_Y":lbm.get_velocity()[:,:,1],
              "Kinetic_Energy":lbm.get_kinetic_energy(),
              }
    
        # output visualization
        if (iter%export_interval==0): 
            visualization.save_pngs(current_time = iter*dt, 
                      export_iter=iter//export_interval, ds=ds, particles=particles, cmap="turbo")
            visualization.save_csvs(export_iter=iter//export_interval, ds=ds)
            
            # postprocessing
            postprocess.run(iter)
        
        # LBM update
        lbm.update(bc_dict)
        if particles.get_n_particles():
            particles.update(lbm.vel)        

    end = time.time()

    ## generate videos from pngs
    visualization.generate_video(ds=ds, fps=ctrl_params["SimSetting"]["fps"], remove_pngs=ctrl_params["SimSetting"]["removePng"])

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