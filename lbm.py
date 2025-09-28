import numpy as np

def stream(width, height, f):
    
    # Stream all internal cells
    n_offset = 1
    for x in range(n_offset, width-1-n_offset):
        for y in range(n_offset, height-1-n_offset):

            f[y*width +x, 1] = f[y*width +x-1, 1]
            f[y*width + x, 2] = f[(y-1)*width + x, 2]
            f[y*width + x-1, 3] = f[y*width + x, 3]
            f[(y-1)*width + x, 4] = f[y*width + x, 4]
            f[(y)*width + (x), 5] = f[(y-1)*width + (x-1), 5]
            f[(y)*width + (x-1), 6] = f[(y-1)*width + (x), 6]
            f[(y-1)*width + (x-1), 7] = f[(y)*width + (x), 7]
            f[(y-1)*width + (x), 8] = f[(y)*width + (x-1), 8]
            
    # Tidy up the edges
    # x = width
    # for y in range(1, height-1):
    #     f[(y)*width + x, 2] = f[(y-1)*width + x, 2]
    #     f[(y-1)*width + x, 4] = f[(y)*width + x, 4]


def collide(width, height, w_, e_, f, rho, u, speed2, wall, omega):
    n_offset = 2
    for x in range(n_offset, width-n_offset):
        for y in range(n_offset, height-n_offset):
            
            i = y*width + x
            
            # Skip over cells containing barriers
            if (wall[i]):
                continue
            else:
                # Compute the macroscopic density
                rho[i] = np.sum(f[i,:])
                # Compute the macroscopic velocities
                if (rho[i] > 0):
                    u[i,:] = np.dot(f[i,1:], e_[1:]) * (1-(rho[i]-1)+((rho[i]-1)**2))
                speed2[i] = np.dot(u[i,:], u[i,:])

                for qi in range(1,9):
                    f_eq = w_[qi] * (1 + 3 * np.dot(u[i,:], e_[qi]) + 4.5 * (np.dot(u[i,:], e_[qi]))**2 - 1.5 * np.dot(u[i,:], u[i,:]))
                    f[i, qi] +=  omega * (f_eq - f[i, qi])
                
                # Conserve mass
                f[i, 0]  = rho[i] - np.sum(f[i,1:])


def bounce(width, height, f, wall):
    n_offset = 2 # original: 2
    # Loop through all interior cells
    for x in range(n_offset, width-n_offset):
        for y in range(n_offset, height-n_offset):
            
            # If the cell contains a boundary
            if (wall[y*width + x]):
                
                # Push densities back 
                f[(y+1)*width+x, 2] = f[y*width+x, 4]
                f[(y-1)*width+x, 4] = f[y*width+x, 2]
                f[y*width+(x-1), 1] = f[y*width+x, 3]
                f[y*width+(x+1), 3] = f[y*width+x, 1]
                f[(y+1)*width+(x+1), 5] = f[y*width+x, 7]
                f[(y+1)*width+(x-1), 6] = f[y*width+x, 8]
                f[(y-1)*width+(x+1), 8] = f[y*width+x, 6]
                f[(y-1)*width+(x-1), 7] = f[y*width+x, 5]
                
                # Clear the densities in the barrier cells
                for fi in f[y*width+x,:]:
                    fi = 0
                
