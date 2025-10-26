from common_modules import np, nb

#@nb.jit(parallel=True, fastmath=True)
def stream(f, e_):
    _, _, q = f.shape

    # streaming
    for i, ei in zip(range(q), e_):
        cx, cy = ei[0], ei[1]
        f[:,1:-1,i] = np.roll(f[:,1:-1,i], cx, axis=1)
        f[:,0,i] = f[:,1,i]
        f[:,-1,i] = f[:,-2,i]
        f[1:-1,:,i] = np.roll(f[1:-1,:,i], cy, axis=0)
        f[0,:,i] = f[1,:,i]
        f[-1,:,i] = f[-2,:,i]

@nb.jit(parallel=True, fastmath=True) 
def bounce(f, wall):
    ny, nx, _ = f.shape

    for yi in range(ny):
        for xi in range(nx):
            
            # If the cell is a wall cell
            if (wall[yi, xi]):
                
                # bounce back
                f[yi+1, xi, 2] = f[yi, xi, 4]
                f[yi-1, xi, 4] = f[yi, xi, 2]
                f[yi, xi+1, 1] = f[yi, xi, 3]
                f[yi, xi-1, 3] = f[yi, xi, 1]
                f[yi+1, xi+1, 5] = f[yi, xi, 7]
                f[yi+1, xi-1, 6] = f[yi, xi, 8]
                f[yi-1, xi+1, 8] = f[yi, xi, 6]
                f[yi-1, xi-1, 7] = f[yi, xi, 5]
                
                # clear f in wall
                f[yi, xi, :] = 0
    