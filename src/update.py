from common_modules import np, nb


#@nb.jit(parallel=True, fastmath=True)
def stream(f, e_):
    ny, nx, q = f.shape
    out = f

    # streaming
    for i, ei in zip(range(q), e_):
        cx, cy = ei[0], ei[1]
        out[:,:,i] = np.roll(out[:,:,i], cx, axis=1)
        out[:,:,i] = np.roll(out[:,:,i], cy, axis=0)

    return out


@nb.jit(parallel=True, fastmath=True) 
def bounce(f, wall):
    ny, nx, _ = f.shape
    f_new = f

    n_offset = 0
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

    