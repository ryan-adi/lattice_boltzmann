from common_modules import np, nb

def stream(f, e_):
    q = len(e_)

    # streaming
    for qi, ei in zip(range(q), e_):
        cx, cy = int(ei[0]), int(ei[1])

        f[:,:,qi] = np.roll(f[:,:,qi], cx, axis=1)
        f[:,:,qi] = np.roll(f[:,:,qi], cy, axis=0)

# # @nb.jit(parallel=True, fastmath=True, cache=True)
# def stream(f, e_):
#     ny, nx, q = f.shape
#     f_old = f.copy()

#     # streaming
#     for x in range(nx):
#             for y in range(ny):
#                 for qi, ei in zip(range(q), e_):
#                     cx, cy = ei[0], ei[1]
#                     x2 = int((x + cx + ny) % nx)
#                     y2 = int((y + cy + ny) % ny)
#                     f[y2, x2, qi] = f_old[y, x, qi] 

@nb.jit(parallel=True, fastmath=True, cache=True)
def bounce(f, wall):
    ny, nx, _ = f.shape

    for yi in nb.prange(ny):
        for xi in nb.prange(nx):
            # If the cell is a wall cell
            if (wall[yi, xi]==1):

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
                f[yi, xi, :] = float('nan')

@nb.jit(parallel=True, fastmath=True, cache=True)
def calc_fe(f, e):
    ny, nx, _ = f.shape
    temp = np.zeros((ny,nx,2))
    for j in nb.prange(ny):
        for i in nb.prange(nx):
            temp[j,i,0] = np.sum(f[j,i,:] * e[:,0])
            temp[j,i,1] = np.sum(f[j,i,:] * e[:,1])
    return temp




   