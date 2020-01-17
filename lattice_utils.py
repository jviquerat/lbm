# Generic imports
import PIL
import copy
import numpy             as np
import matplotlib        as mplt
import matplotlib.pyplot as plt

### ************************************************
### Class defining lattice object
class Lattice:
    ### ************************************************
    ### Constructor
    def __init__(self, name,
                 xmin, xmax,
                 ymin, ymax,
                 nx,   ny,
                 q,    tau):

        self.name = name
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.nx   = nx
        self.ny   = ny
        self.q    = q
        self.tau  = tau

    ### ************************************************
    ### Solve LBM
    def solve(self, it_max, rho):

        self.g    = np.zeros((self.q,  self.nx, self.ny))
        self.g_eq = np.zeros((self.q,  self.nx, self.ny))
        self.g_up = np.zeros((self.q,  self.nx, self.ny))
        self.rho  = np.ones ((self.nx, self.ny))*rho
        self.u    = np.zeros((2,       self.nx, self.ny))

        # Missing initial velocity profile here

        # Solve equilibrium to get a starting distribution
        for q in range(self.q):
                self.g[q,:,:] = self.rho[:,:]*self.w[q]*(
                    1.0 + 3.0*(self.u[0,:,:]*self.c[q,0] +
                               self.u[1,:,:]*self.c[q,1]))

        # Solve
        for it in range(it_max):

            # Streaming
            for q in range(self.q):
                self.g_up = np.roll(np.roll(
                    self.g[q,:,:],self.c[q,0],axis=0),
                    self.c[q,1],axis=1)

            # Update rho
            self.rho = np.sum(self.g_up, axis=0)

            # Update u
            self.u[:,:,:] = 0.0

            for q in range(self.q):
                self.u[0,:,:] += self.c[q,0]*self.g_up[q,:,:]
                self.u[1,:,:] += self.c[q,1]*self.g_up[q,:,:]

            self.u[0,:,:] /= self.rho[:,:]
            self.u[1,:,:] /= self.rho[:,:]

            # Equilibrium
            for q in range(self.q):
                self.g_eq[q,:,:] = self.rho[:,:]*self.w[q]*(
                    1.0 + 3.0*(self.u[0,:,:]*self.c[q,0] +
                               self.u[1,:,:]*self.c[q,1]))

            # Collisions
            self.g = self.g_up - (1.0/self.tau)*(self.g_up - self.g_eq)

            # Inflow  b.c.


            # Outflow b.c.
            self.g[right,-1,:] = self.g[right,-2,:]

            # Top/bottom walls b.c.

            # Obstacle b.c.


    ### ************************************************
    ### Solve equilibrium state
    #def solve_eq(self, rho, u):

    ### ************************************************
    ### Initialize computation
    def init_computation(self, u_in):

        # D2Q9 Velocities
        self.c = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]])
        print(self.c)

        # Weights
        # Cardinal values, then extra-cardinal values, then central value
        idx_card       = [np.linalg.norm(ci)<1.1 for ci in self.c]
        idx_extra_card = [np.linalg.norm(ci)>1.1 for ci in self.c]

        self.w                             = np.ones(self.q)
        self.w[np.asarray(idx_card)]       = 1./9.
        self.w[np.asarray(idx_extra_card)] = 1./36.
        self.w[0]                          = 4./9.

        # Boundary conditions
        # Velocities on which to apply the different BC
        self.right = np.arange(self.q)[np.asarray([ci[0] <0
                                                   for ci in self.c])]
        self.left  = np.arange(self.q)[np.asarray([ci[0] >0
                                                   for ci in self.c])]
        self.mid   = np.arange(self.q)[np.asarray([ci[0]==0
                                                   for ci in self.c])]
        #self.ns    = [self.c.tolist().index((-self.c[i]).tolist()) for i in range(self.q)]
        self.ns    = self.c.copy()

        # Initial velocity
        #self.u_in  = np.fromfunction(lambda d,x,y: (1-d)*uLB,(2,nx,ny))

    ### ************************************************
    ### Generate lattice
    def generate(self, polygon):

        # Because we loop on the lattice left-right and top-down,
        # we need to flip the polygon up-down
        poly       = polygon.copy()
        poly[:,1] *= -1.0

        # Declare lattice array
        self.lattice = np.zeros((self.ny, self.nx), dtype=bool)

        # Fill lattice
        for i in range(self.nx):
            for j in range(self.ny):
                pt           = self.lattice_coords(j, i)
                inside       = self.is_inside(poly, pt)
                self.lattice[j,i] = inside

    ### ************************************************
    ### Get lattice coordinates from integers
    def lattice_coords(self, j, i):

        # Compute and return the coordinates of the lattice node (i,j)
        dx = (self.xmax - self.xmin)/(self.nx - 1)
        dy = (self.ymax - self.ymin)/(self.ny - 1)
        x  = self.xmin + i*dx
        y  = self.ymin + j*dy

        return [x, y]

    ### ************************************************
    ### Determine if a pt is inside or outside a closed polygon
    def is_inside(self, poly, pt):

        # Initialize
        j         = len(poly) - 1
        odd_nodes = False

        # Check if point is inside or outside
        # This is a valid algorithm for any non-convex polygon
        for i in range(len(poly)):
            if (((poly[i,1] < pt[1] and poly[j,1] >= pt[1])  or
                 (poly[j,1] < pt[1] and poly[i,1] >= pt[1])) and
                 (poly[i,0] < pt[0] or  poly[j,0]  < pt[0])):

                # Compute slope
                slope = (poly[j,0] - poly[i,0])/(poly[j,1] - poly[i,1])

                # Check side
                if ((poly[i,0] + (pt[1] - poly[i,1])*slope) < pt[0]):
                    odd_nodes = not odd_nodes

            # Increment
            j = i

        return odd_nodes

    ### ************************************************
    ### Generate lattice image
    def generate_image(self):

        # Plot and save image of lattice
        plt.axis('off')
        plt.imshow(self.lattice,
                   cmap=mplt.cm.inferno)
        plt.savefig(self.name, dpi=200, bbox_inches='tight')
        plt.close()
        self.trim_white(self.name+'.png')

    ### ************************************************
    ### Crop white background from image
    def trim_white(self, filename):

        # Trim using PIL
        im   = PIL.Image.open(filename)
        bg   = PIL.Image.new(im.mode, im.size, (255,255,255))
        diff = PIL.ImageChops.difference(im, bg)
        bbox = diff.getbbox()
        cp   = im.crop(bbox)
        cp.save(filename)
