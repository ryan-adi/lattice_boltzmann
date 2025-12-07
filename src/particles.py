from common_modules import np

def isInside(x, xmin, xmax):
        return (xmin+1<x) & (x<xmax-1)

class Particles():
    '''Defines class for particles that move w.r.t. flow field but doesn't interact with each other.'''
    def __init__(self, x_bound:list, y_bound:list):
        self.position = np.array([])
        self.xmin = x_bound[0]
        self.ymin = y_bound[0]
        self.xmax = x_bound[1]
        self.ymax = y_bound[1]

    def initialize(self, params:dict):
        for key, val in params.items():
            if key=="spawnGrid":
                for i in range(val.shape[0]):
                    x_seq = val[i][:3]
                    y_seq = val[i][3:]
                    self.spawn_grid(x_seq, y_seq)
            elif key=="spawnCircular":
                for i in range(val.shape[0]):
                    n_particle = val[i][0]
                    radius = val[i][1]
                    center = val[i][2:]
                    self.spawn_circular(n_particle, radius, center)
            else:
                Warning("Spawning method not found.")

    def spawn_grid(self, x_seq:list, y_seq:list):
        y_1d = np.arange(y_seq[0], y_seq[1], y_seq[2])
        x_1d = np.arange(x_seq[0], x_seq[1], x_seq[2])
        xy = np.array([[x,y] for x in x_1d for y in y_1d], dtype=float)
        if self.position.size==0:
            self.position = xy
        else:
            self.position = np.concatenate((self.position, xy), axis=0)

    def spawn_circular(self, n_particle, radius, center:list):
        theta = np.linspace(0, 359.999, n_particle) * np.pi / 2
        x_coord = radius * np.cos(theta) + center[0]
        y_coord = radius * np.sin(theta) + center[1]
        xy = np.stack([x_coord, y_coord]).T
        if self.position.size==0:
            self.position = xy
        else:
            self.position = np.concatenate((self.position, xy), axis=0)

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
        return self.position.shape[0]