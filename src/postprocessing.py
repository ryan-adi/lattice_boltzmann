from common_modules import plt, os, subprocess, glob, shutil, np

def set_colorbar(im, label):
    cbar = plt.colorbar(im)
    plt.gca().invert_yaxis()
    cbar.set_label(label, rotation=270)

def create_output_folder(case_name):
    case_pngs = os.path.join("output", case_name, "pngs")
    case_csvs = os.path.join("output", case_name, "csvs")

    if os.path.exists(case_pngs):
        # if already exist delete case dir
        case_dir = os.path.join("output", case_name)
        shutil.rmtree(case_dir)
        
    os.makedirs(case_pngs)
    os.makedirs(case_csvs)

def save_csvs(case_name, current_iter, ds:dict):
     
    # define main directory
    main_dir = os.getcwd()
    case_csvs = os.path.join("output", case_name, "csvs")
    
    os.chdir(case_csvs)
    # create csvs for all variables in ds
    for var, data in ds.items():

        var_name = ''.join([c for c in var if c.isupper()])
        csv_name = f"t{current_iter:06d}_{var_name}.csv" 

        if not os.path.exists(var_name):
            os.mkdir(var_name)
        os.chdir(var_name)
        
        # filter out zero entries
        data[data < 1e-6] = 0

        np.savetxt(csv_name, data, delimiter=",", fmt="%.3e")

        os.chdir("..")



    # go to directory where main is
    os.chdir(main_dir)

def save_pngs(case_name, current_iter, ds:dict):
    
    # define main directory
    main_dir = os.getcwd()
    case_pngs = os.path.join("output", case_name, "pngs")
    
    os.chdir(case_pngs)
    # create pngs for all variables in ds
    for var, data in ds.items():
        fig, ax = plt.subplots(figsize=(20,5))

        im = plt.imshow(data, vmin=0.0, vmax=0.05, cmap='viridis') #  
        set_colorbar(im, label=var)

        var_name = ''.join([c for c in var if c.isupper()])
        png_name = f"t{current_iter:06d}_{var_name}.png" 

        if not os.path.exists(var_name):
            os.mkdir(var_name)
        os.chdir(var_name)

        fig.savefig(png_name)
        
        os.chdir("..")

        plt.close()

    # go to directory where main is
    os.chdir(main_dir)

def generate_video(case_name:str, ds:dict, fps=20, remove_pngs=True):
    # define pngs folder
    main_dir = os.getcwd()
    case_pngs = os.path.join("output", case_name, "pngs")
    os.chdir(case_pngs)

    # call ffmpeg to convert pngs to mp4
    for var, _ in ds.items():
        # file names
        var_name = ''.join([c for c in var if c.isupper()])
        video_name= case_name + "_" + var_name + ".mp4"
        png_name = "t%06d_" + var_name + ".png"

        os.chdir(var_name)

        # call ffmpeg
        fps = str(fps)
        subprocess.call([
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