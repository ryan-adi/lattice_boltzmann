from common_modules import plt, os, subprocess, glob, shutil, np
from src.lbm import * 

class Visualization():
    def __init__(self, lbm:LatticeBoltzmann, config:dict):
        self.case_name = config["CaseName"]
        
        # plot settings
        plt.rcParams["figure.figsize"] = (config["figSize"][0],config["figSize"][1])
        plt.rcParams.update({
            'axes.labelsize': 15,      # Font size for axis labels
            'axes.titlesize': 16,      # Font size for axis titles
            'xtick.labelsize': 12,     # Font size for x-axis ticks
            'ytick.labelsize': 12,     # Font size for y-axis ticks
            'font.family': 'serif',    # Font family (default is 'sans-serif')
            'font.size': 12,           # General font size for text in the plot
        })

    def set_colorbar(self, im, label):
        cbar = plt.colorbar(im)
        plt.gca().invert_yaxis()
        cbar.set_label(label, fontweight='bold', labelpad=35, rotation=270)

    def create_output_folder(self):
        case_pngs = os.path.join("output", self.case_name, "visuals")
        case_csvs = os.path.join("output", self.case_name, "csvs")
            
        os.makedirs(case_pngs)
        os.makedirs(case_csvs)

    def axis_limits(self, var):
        match var:
            case "Density":
                value_min = 0.0
                value_max = 2.1
            case "Velocity_X":
                value_min = -0.3
                value_max = 0.3
            case "Velocity_Y":
                value_min = -0.3
                value_max = 0.3
            case "Kinetic_Energy":
                value_min = 0.0
                value_max = 0.1
            case "Scalar":
                value_min = 0.0
                value_max = 5.0
            case _:
                value_min = None
                value_max = None

        return value_min, value_max

    def save_csvs(self, export_iter, ds:dict):
        
        # define main directory
        main_dir = os.getcwd()
        case_csvs = os.path.join("output", self.case_name, "csvs")
        
        os.chdir(case_csvs)
        # create csvs for all variables in ds
        for var, data in ds.items():

            var_name = var #''.join([c for c in var if c.isupper()])
            csv_name = f"t{export_iter:06d}_{var_name}.csv" 

            if not os.path.exists(var_name):
                os.mkdir(var_name)
            os.chdir(var_name)
            
            # flip data to in y-axis
            data = np.flip(data,0) 

            # filter out zero entries
            # if (data.min() >= 0.0):
            #     data[data < 1e-6] = 0.0

            np.savetxt(csv_name, data, delimiter=",", fmt="%.3e")

            os.chdir("..")

        # go to directory where main is
        os.chdir(main_dir)

    def save_pngs(self, current_time, export_iter, ds:dict, particles, cmap):
        
        # define main directory
        main_dir = os.getcwd()
        case_pngs = os.path.join("output", self.case_name, "visuals")
        os.chdir(case_pngs)

        # create pngs for all variables in ds
        for var, data in ds.items():
            fig, ax = plt.subplots()
            ax.set_xlabel("x Coordinates")
            ax.set_ylabel("y Coordinates")
            timestamp = f"t={current_time:.1f} s"
            ax.text(0.0, 1.1, timestamp, 
                    fontsize=23, fontweight='bold',
                    bbox=dict(facecolor='red', alpha=0.5), 
                    transform=ax.transAxes)

            value_min = value_max = 0.0
            value_min, value_max = self.axis_limits(var)

            im = plt.imshow(data,  cmap=cmap)
            im = plt.imshow(data,  vmin=value_min, vmax=value_max, cmap=cmap)
            self.set_colorbar(im, label=var)

            var_name = var
            png_name = f"t{export_iter:06d}_{var_name}.png" 

            # create points in plot
            particle_pos = particles.get_position()
            if particle_pos.size:
                plt.scatter(particle_pos[:,0], particle_pos[:,1], c="red", edgecolors="k")

            # save png
            if not os.path.exists(var_name):
                os.mkdir(var_name)
            os.chdir(var_name)
            fig.savefig(png_name)
            os.chdir("..")

            plt.close()

        # go to directory where main is
        os.chdir(main_dir)

    def generate_video(self, ds:dict, fps:int, remove_pngs=True):
        # define pngs folder
        main_dir = os.getcwd()
        case_pngs = os.path.join("output", self.case_name, "visuals")
        os.chdir(case_pngs)

        # call ffmpeg to convert pngs to mp4
        for var, _ in ds.items():
            # file names
            var_name = var #''.join([c for c in var if c.isupper()])
            video_name= self.case_name + "_" + var_name + ".mp4"
            png_name = "t%06d_"+ var_name +".png" 

            os.chdir(var_name)

            # call ffmpeg
            fps = str(fps)
            subprocess.run([
                'ffmpeg.exe', '-framerate', fps, 
                '-i', png_name, 
                '-r', '30', 
                '-pix_fmt', 'yuv420p',
                video_name])
            
            # move mp4 to case dir
            shutil.move(video_name, "../"+video_name)

            # remove pngs
            if remove_pngs:
                for file_name in glob.glob("*_" + var_name + ".png"):
                    os.remove(file_name)

            os.chdir("..")

        # go to directory where main is
        os.chdir(main_dir)