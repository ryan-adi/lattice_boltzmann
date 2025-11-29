from common_modules import np, plt, os, shutil

# plot settings
plt.rcParams["figure.figsize"] = (50,3)
plt.rcParams.update({
    'axes.labelsize': 14,      # Font size for axis labels
    'axes.titlesize': 16,      # Font size for axis titles
    'xtick.labelsize': 12,     # Font size for x-axis ticks
    'ytick.labelsize': 12,     # Font size for y-axis ticks
    #'font.family': 'serif',    # Font family (default is 'sans-serif')
    'font.size': 12,           # General font size for text in the plot
    'lines.linewidth':5         # line width
})

class Postprocess():
    def __init__(self, lb, case_name):
        self.lb = lb
        self.case_name = case_name
        self.iter = 0

    def create_postprocessing_folder(self):
        case_postprocessing = os.path.join("output", self.case_name, "postprocessing")      
        os.makedirs(case_postprocessing)

    def plot_cross_section_value(self, variable, loc):
        var_name = "Velocity"

        main_dir = os.getcwd()
        case_postprocessing = os.path.join("output", self.case_name, "postprocessing")
        os.chdir(case_postprocessing)

        x = np.arange(0, self.lb.ny, 1)
        y = variable[:, loc]

        fig, ax = plt.subplots(figsize=(20,5))
        ax.set_xlabel("x Coordinates")
        ax.set_ylabel(var_name)
        ax.set_xlim(0,self.lb.ny)
        ax.set_ylim(0,0.3)
        x_mean = np.sum(x*y)/np.sum(y)
        ax.vlines(x_mean, 0.0, 0.5, "r")
        # ax.annotate(str(x_mean), xy =(x_mean, 0.25),
        #      xytext =(x_mean+1, 0.25),
        #      arrowprops = dict(facecolor ='green',
        #                        shrink = 0.05),   )
        # timestamp = f"t={self.time:.1f} s"
        # ax.text(0.0, 1.1, timestamp, 
        #         fontsize=23, fontweight='bold',
        #         bbox=dict(facecolor='red', alpha=0.5), 
        #         transform=ax.transAxes)
        plt.plot(x, y)

        png_name = f"slice_{var_name}_{self.iter}.png" 
        fig.savefig(png_name)
        plt.close()
        
        os.chdir(main_dir)

    def run(self, iter):
        self.iter = iter
        self.plot_cross_section_value(self.lb.u[:,:,0], loc=2)