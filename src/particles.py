from common_modules import np

def isInside(x, xmin, xmax):
        return (xmin<x) & (x<xmax)

class Particles:
    '''Defines class for particles that move w.r.t. flow field but doesn't interact with each other.'''
    def __init__(self, pos_0:np.array, xmax, ymax):
        self.position = pos_0.astype(float)
        self.n_p = pos_0.shape[0]
        self.xmin = 1
        self.ymin = 1
        self.xmax = xmax-1
        self.ymax = ymax-1

    def calc_position(self, u_field):
        for pi in range(self.position.shape[0]):
            x,y = int(self.position[pi][0]), int(self.position[pi][1])
            stencil = np.array([[x,y], [x+1,y], [x,y+1], [x-1,y], [x,y-1]])
            vel = 0.0
            for id_i in stencil:
                vel += u_field[id_i[1], id_i[0]]

            self.position[pi,0] += vel[0] / stencil.shape[0]
            self.position[pi,1] += vel[1] / stencil.shape[0]
        
    def despawn(self):
        temp = self.position.copy()
        mask = isInside(temp[:,0], self.xmin, self.xmax) & isInside(temp[:,1], self.ymin, self.ymax)
        temp = temp[mask]
        self.position = temp
      

    def update(self, u_field):
        self.calc_position(u_field)
        self.despawn()
        
    # ========== GET FUNCTIONS ========== #
    def get_position(self):
        return self.position
    
    def get_n_particles(self):
        return self.n_p