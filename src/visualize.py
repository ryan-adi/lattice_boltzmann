from common_modules import plt, os, subprocess, glob, shutil, imageio

def set_colorbar(im, label):
    cbar = plt.colorbar(im)
    plt.gca().invert_yaxis()
    cbar.set_label(label, rotation=270)

def create_output_folder(case_name):
    case_pngs = os.path.join("output", case_name, "pngs")
    if not os.path.exists(case_pngs):
        os.makedirs(case_pngs)
    else:
        case_dir = os.path.join("output", case_name)
        shutil.rmtree(case_dir)
        os.makedirs(case_pngs)

def save_pngs(case_name, current_iter, ds:dict):
    
    # define main directory
    main_dir = os.getcwd()
    case_pngs = os.path.join("output", case_name, "pngs")
    
    os.chdir(case_pngs)
    # create pngs for all variables in ds
    fig, ax = plt.subplots(figsize=(20,5))
    for var, data in ds.items():

        im = plt.imshow(data, vmin=0.0, vmax=0.05, cmap='viridis') #  
        set_colorbar(im, label=var)

        var_name = ''.join([c for c in var if c.isupper()])
        png_name = f"t{current_iter:06d}_{var_name}.png" 
        fig.savefig(png_name)
        plt.close()

    # go to directory where main is
    os.chdir(main_dir)

def generate_video(case_name:str, ds:dict, fps=20, save_pngs=True):
    # define pngs folder
    main_dir = os.getcwd()
    case_pngs = os.path.join("output", case_name, "pngs")
    os.chdir(case_pngs)

    # call ffmpeg to convert pngs to mp4
    for var, _ in ds.items():
        # 
        var_name = ''.join([c for c in var if c.isupper()])
        video_name= case_name + "_" + var_name + ".mp4"
        png_name = "t%06d_" + var_name + ".png"

        # call ffmpeg
        fps = str(fps)
        subprocess.call([
            'ffmpeg.exe', '-framerate', fps, 
            '-i', png_name, 
            '-r', '30', 
            '-pix_fmt', 'yuv420p',
            video_name])

        # remove pngs
        if save_pngs:
            for file_name in glob.glob("*_" + var_name + ".png"):
                os.remove(file_name)

    # go to directory where main is
    os.chdir(main_dir)