from common_modules import np, nb


#@nb.jit
def stream(f):
    ny, nx, q = f.shape
    out = f
        
    # n_offset = 1
    # for yi in range(n_offset, ny-n_offset):
    #         for xi in range(n_offset, nx-n_offset):
    #             for k in range(1,q):
    #                 if (k<5):
    #                     direction = 0.5*(k-1)*np.pi
    #                 else:
    #                     direction = (0.5*(k-5)+0.25)*np.pi

    #                 x0 = int(np.round(xi-np.cos(direction)))
    #                 y0 = int(np.round(yi-np.sin(direction)))

    #                 out[yi,xi,k] = f[y0, x0 ,k]

    # out[:,1:nx,1] = f[:,0:nx-1,1]
    # out[1:ny,:,2] = f[0:ny-1,:,2]
    # out[:,0:nx-1,3] = f[:,1:nx,3]
    # out[0:ny-1,:,4] = f[1:ny,:,4]
    # out[1:ny,1:nx,5] = f[0:ny-1,0:nx-1,5]
    # out[1:ny,0:nx-1,6] = f[0:ny-1,1:nx,6]
    # out[0:ny-1,0:nx-1,7] = f[1:ny,1:nx,7]
    # out[0:ny-1,1:nx,8] = f[1:ny,0:nx-1,8]

    e_ = np.array([np.array([0.0, 0.0]),
                            np.array([1.0, 0.0]),np.array([0.0, 1.0]),
                            np.array([-1.0, 0.0]),np.array([0.0, -1.0]),
                            np.array([1.0, 1.0]),np.array([-1.0, 1.0]),
                            np.array([-1.0, -1.0]),np.array([1.0, -1.0])])
    # streaming
    for i, ei in zip(range(9), e_):
        cx, cy = ei[0], ei[1]
        out[:,:,i] = np.roll(out[:,:,i], cx, axis=1)
        out[:,:,i] = np.roll(out[:,:,i], cy, axis=0)

    return out


#@nb.jit 
def bounce(f, wall):
    ny, nx, _ = f.shape
    f_new = f

    n_offset = 2 
    for yi in range(n_offset, ny-n_offset):
        for xi in range(n_offset, nx-n_offset):
            
            # If the cell is a wall cell
            if (wall[yi, xi]):
                
                # bounce back
                f_new[yi+1, xi, 2] = f[yi, xi, 4]
                f_new[yi-1, xi, 4] = f[yi, xi, 2]
                f_new[yi, xi+1, 1] = f[yi, xi, 3]
                f_new[yi, xi-1, 3] = f[yi, xi, 1]
                f_new[yi+1, xi+1, 5] = f[yi, xi, 7]
                f_new[yi+1, xi-1, 6] = f[yi, xi, 8]
                f_new[yi-1, xi+1, 8] = f[yi, xi, 6]
                f_new[yi-1, xi-1, 7] = f[yi, xi, 5]
                
                # clear f in wall
                f_new[yi, xi, :] = 0

    return f_new

    