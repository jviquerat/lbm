# Generic imports
import os
import PIL
import progress.bar
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
                 q,    tau,
                 output_dir):

        self.name       = name
        self.xmin       = xmin
        self.xmax       = xmax
        self.ymin       = ymin
        self.ymax       = ymax
        self.nx         = nx
        self.ny         = ny
        self.q          = q
        self.tau        = tau
        self.output_dir = output_dir
        self.png_dir    = self.output_dir+'./png/'
        self.output_it  = 0

        if (not os.path.exists(self.output_dir)): os.makedirs(self.output_dir)
        if (not os.path.exists(self.png_dir)):    os.makedirs(self.png_dir)

    ### ************************************************
    ### Solve LBM
    def solve(self, it_max, u_in, rho, freq):

        self.g    = np.zeros((self.q,  self.ny, self.nx))
        self.g_eq = np.zeros((self.q,  self.ny, self.nx))
        self.g_up = np.zeros((self.q,  self.ny, self.nx))
        self.rho  = np.ones ((self.ny, self.nx))*rho
        self.u    = np.zeros((2,       self.ny, self.nx))

        # Input velocity profile
        self.u_in  = u_in*np.fromfunction(self.poiseuille,
                                          (2, self.ny, self.nx))

        #self.u[:,:,0] = self.u_in[:,:,0]
        self.u                 = self.u_in
        self.u[:,self.lattice] = 0.0

        # Initial distribution
        self.equilibrium(self.g, self.rho, self.u)

        # Solve with progress bar
        bar = progress.bar.Bar('Solving...', max=it_max)
        for it in range(it_max+1):

            # Sizes
            lx = self.nx-1
            ly = self.ny-1

            # Compute macroscopic fields ### local
            self.macro()

            # Inflow b.c. : Zou-He, macro part ### local
            self.u[:,:,0] = self.u_in[:,:,0]
            self.rho[:,0] = 1.0/(1.0-self.u[0,:,0])*( \
                np.sum(self.g[self.mid, :,0],axis=0)  \
                + 2.0*np.sum(self.g[self.left,:,0],axis=0))

            # Outflow b.c. : Zou-He, macro part ### local
            self.u[1,:,lx] = 0.0
            self.rho[:,lx] = 1.0
            self.u[0,:,lx] =-1.0 + np.sum(self.g[self.mid,:,lx],axis=0) \
                + 2.0*np.sum(self.g[self.right,:,lx],axis=0)

            # Obstacle b.c. : macro part ### local but need to cut and distribute lattice
            self.u[:,self.lattice] = 0.0

            # Compute equilibrium state ### local
            self.equilibrium(self.g_eq, self.rho, self.u)

            # Collisions ### local
            self.g_up = self.g - (1.0/self.tau)*(self.g - self.g_eq)

            # Streaming ### NOT local : use extended arrays for artificial boundaries
                        ###             and then just stream normally ?

            # Stream 0
            self.g[0, :, :]           = self.g_up[0, :, :]

            # Stream +x and -x
            self.g[1, 1:ly-1, 1:lx]   = self.g_up[1, 1:ly-1, 0:lx-1]
            self.g[2, 1:ly-1, 0:lx-1] = self.g_up[2, 1:ly-1, 1:lx]

            # Stream +y and -y
            self.g[3, 0:ly-1, 1:lx-1] = self.g_up[3, 1:ly,   1:lx-1]
            self.g[4, 1:ly,   1:lx-1] = self.g_up[4, 0:ly-1, 1:lx-1]

            # Stream +x+y and -x-y
            self.g[5, 0:ly-1, 1:lx]   = self.g_up[5, 1:ly,   0:lx-1]
            self.g[6, 1:ly,   0:lx-1] = self.g_up[6, 0:ly-1, 1:lx]

            # Stream -x+y and +x-y
            self.g[7, 0:ly-1, 0:lx-1] = self.g_up[7, 1:ly,   1:lx]
            self.g[8, 1:ly,   1:lx]   = self.g_up[8, 0:ly-1, 0:lx-1]

            # Inflow b.c. : Zou-He, micro part ### local
            self.g[1,:,0] = self.g_eq[1,:,0] + self.g[2,:,0] - self.g_eq[2,:,0]
            self.g[5,:,0] = 0.5*(self.rho[:,0]*self.u[0,:,0] - self.g[1,:,0] \
                + self.g[2,:,0] - self.g[3,:,0] + self.g[4,:,0] + 2.0*self.g[6,:,0])
            self.g[8,:,0] = self.g[3,:,0] - self.g[4,:,0] + self.g[5,:,0] \
                - self.g[6,:,0] + self.g[7,:,0]

            # Outflow b.c. : Zou-He, micro part ### local
            self.g[2,:,lx] = self.g_eq[2,:,lx] + self.g[1,:,lx] - self.g_eq[1,:,lx]
            self.g[6,:,lx] =-0.5*(self.rho[:,lx]*self.u[0,:,lx] - self.g[1,:,lx] \
                + self.g[2,:,lx] - self.g[3,:,lx] + self.g[4,:,lx] - 2.0*self.g[5,:,lx])
            self.g[7,:,lx] =-self.g[3,:,lx] + self.g[4,:,lx] - self.g[5,:,lx] \
                + self.g[6,:,lx] + self.g[8,:,lx]

            # Top b.c. : Zou-He, micro part ### local
            self.g[4,0,:] = self.g_eq[4,0,:] + self.g[3,0,:] - self.g_eq[3,0,:]
            self.g[6,0,:] = 0.5*(self.g[1,0,:] + self.g[3,0,:] - self.g[2,0,:] \
                - self.g[4,0,:] + 2.0*self.g[5,0,:])
            self.g[8,0,:] = self.g[3,0,:] - self.g[4,0,:] + self.g[5,0,:] \
                - self.g[6,0,:] + self.g[7,0,:]

            # Bottom b.c. : Zou-He, micro part ### local
            self.g[3,ly,:] = self.g_eq[3,ly,:] + self.g[4,ly,:] - self.g_eq[4,ly,:]
            self.g[5,ly,:] = 0.5*(2.0*self.g[6,ly,:] - self.g[1,ly,:] - self.g[3,ly,:] \
                + self.g[2,ly,:] + self.g[4,ly,:])
            self.g[7,ly,:] = self.g[8,ly,:] - self.g[3,ly,:] + self.g[4,ly,:] \
                - self.g[5,ly,:] + self.g[6,ly,:]

            # Obstacle b.c. ### local
            for q in range(self.q):
                self.g[q,self.lattice] = self.g[self.ns[q], self.lattice]

            # Output view ### NOT local : need to reconstruct before output
            if (it%freq==0): # Visualization
                plt.clf()
                plt.imshow(np.sqrt(self.u[0]**2+self.u[1]**2),
                           cmap = 'Blues',
                           vmin = 0.0,
                           vmax = u_in)
                plt.colorbar()
                plt.savefig(self.png_dir+"vel."+str(self.output_it)+".png")
                self.output_it += 1

            # Increment bar
            bar.next()

        # End bar
        bar.finish()

    ### ************************************************
    ### Compute equilibrium state
    def equilibrium(self, g, rho, u):

        # Reset distribution
        g[:,:,:] = 0.0

        # Compute velocity term
        v = (3.0/2.0)*(u[0,:,:]**2 + u[1,:,:]**2)

        # Compute equilibrium
        for q in range(self.q):
            t        = 3.0*(u[0,:,:]*self.c[q,0] + u[1,:,:]*self.c[q,1])
            g[q,:,:] = rho*self.w[q]*(1.0 + t + 0.5*t**2 - v)

    ### ************************************************
    ### Compute macroscopic fields
    def macro(self):

        # Compute density
        self.rho = np.sum(self.g, axis=0)

        # Compute velocity
        self.u[:,:,:] = 0.0

        for q in range(self.q):
            self.u[0,:,:] += self.c[q,0]*self.g[q,:,:]
            self.u[1,:,:] += self.c[q,1]*self.g[q,:,:]

        self.u[0,:,:] /= self.rho[:,:]
        self.u[1,:,:] /= self.rho[:,:]

    ### ************************************************
    ### Initialize computation
    def init_computation(self):

        # D2Q9 Velocities
        self.c  = np.array([ [ 0, 0],
                             [ 1, 0], [-1, 0],
                             [ 0, 1], [ 0,-1],
                             [ 1, 1], [-1,-1],
                             [-1, 1], [ 1,-1]])

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
        self.right = np.arange(self.q)[np.asarray([ci[0] >0
                                                   for ci in self.c])]
        self.left  = np.arange(self.q)[np.asarray([ci[0] <0
                                                   for ci in self.c])]
        self.mid   = np.arange(self.q)[np.asarray([ci[0]==0
                                                   for ci in self.c])]
        self.top   = np.arange(self.q)[np.asarray([ci[1] >0
                                                   for ci in self.c])]
        self.bot   = np.arange(self.q)[np.asarray([ci[1] <0
                                                   for ci in self.c])]
        self.ns    = np.asarray([self.c.tolist().index(
            (-self.c[i]).tolist()) for i in range(self.q)])

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
        bar = progress.bar.Bar('Generating...', max=self.nx*self.ny)
        for i in range(self.nx):
            for j in range(self.ny):
                pt           = self.lattice_coords(j, i)
                inside       = self.is_inside(poly, pt)
                self.lattice[j,i] = inside

                bar.next()
        bar.finish()

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

    ### ************************************************
    ### Poiseuille flow
    def poiseuille(self, d, y, x):

        return (1.0-d)*(4.0*y/self.ny)*(1.0 - y/self.ny)
