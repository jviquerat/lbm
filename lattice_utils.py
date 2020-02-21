# Generic imports
import os
import progress.bar
import numpy             as np
import matplotlib        as mplt
import matplotlib.pyplot as plt
import PIL
from   PIL               import Image

### ************************************************
### Class defining lattice object
class Lattice:
    ### ************************************************
    ### Constructor
    def __init__(self,       name,
                 xmin,       xmax,
                 ymin,       ymax,
                 nx,         ny,
                 q,          tau,
                 output_dir, dpi):

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
        self.dpi        = dpi

        if (not os.path.exists(self.output_dir)): os.makedirs(self.output_dir)
        if (not os.path.exists(self.png_dir)):    os.makedirs(self.png_dir)

    ### ************************************************
    ### Initialize lattice
    def init_lattice(self, rho):

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
        c          = self.c
        self.right = np.arange(self.q)[np.asarray([ci[0] >0 for ci in c])]
        self.left  = np.arange(self.q)[np.asarray([ci[0] <0 for ci in c])]
        self.mid   = np.arange(self.q)[np.asarray([ci[0]==0 for ci in c])]
        self.top   = np.arange(self.q)[np.asarray([ci[1] >0 for ci in c])]
        self.bot   = np.arange(self.q)[np.asarray([ci[1] <0 for ci in c])]
        self.ns    = np.asarray([self.c.tolist().index(
            (-self.c[i]).tolist()) for i in range(self.q)])

        # Allocate arrays
        self.g    = np.zeros((self.q,  self.ny, self.nx))
        self.g_eq = np.zeros((self.q,  self.ny, self.nx))
        self.g_up = np.zeros((self.q,  self.ny, self.nx))
        self.rho  = np.ones ((self.ny, self.nx))*rho
        self.u    = np.zeros((2,       self.ny, self.nx))

    ### ************************************************
    ### Solve LBM
    def solve(self, it_max, u_in, rho,
              U_ref, L_ref, dx, dt, freq):

        # Initialize lattice
        self.init_lattice(rho)

        # Input velocity profile
        self.u_in = u_in*np.fromfunction(self.poiseuille,(2,self.ny,self.nx))
        self.u                 = self.u_in
        self.u[:,self.lattice] = 0.0

        # Initial distribution
        self.equilibrium(self.g, self.rho, self.u)
        self.macro()

        # Solve
        bar = progress.bar.Bar('Solving...', max=it_max)
        for it in range(it_max+1):

            # Drag and lift
            self.drag_lift(it, rho, U_ref, L_ref, dx, dt)

            # Compute equilibrium state
            self.equilibrium(self.g_eq, self.rho, self.u)

            # Collisions
            self.g_up = self.g - (1.0/self.tau)*(self.g - self.g_eq)

            # Streaming
            self.stream()

            # Boundary conditions
            self.zou_he_inlet(self.g, self.g_eq, self.rho, self.u)
            self.zou_he_outlet(self.g, self.g_eq, self.rho, self.u)
            self.zou_he_top_wall(self.g)
            self.zou_he_bottom_wall(self.g)
            self.bounce_back_obstacle(self.g)

            # Compute macroscopic fields
            self.macro()

            # Output view
            self.output_view(it, freq, u_in)

            # Increment bar
            bar.next()

        # End bar
        bar.finish()

    ### ************************************************
    ### Compute drag and lift
    def drag_lift(self, it, rho, U_ref, L_ref, dx, dt):

        # Initialize
        force = np.zeros((2))

        # Loop over obstacle array
        for k in range(len(self.obstacle)):
            i = self.obstacle[k,0]
            j = self.obstacle[k,1]

            for q in range(1,self.q):
                qb        = self.ns[q]
                dc        = self.c[q, :]
                ic        = self.c[qb,:]
                ii        = i + dc[0]
                jj        = j + dc[1]
                w         = self.lattice[jj,ii]
                df        = self.g[q,j,i] - self.g[qb,jj,ii]
                force[:] += dc[:]*(1.0-w)*df

        # Normalize coefficient
        force *= dx/dt
        force *= 1.0/(rho*L_ref*U_ref**2)

        # Write to file
        filename = self.output_dir+'drag_lift'
        with open(filename, 'a') as f:
            f.write('{} {} {}\n'.format(it, force[0], force[1]))

    ### ************************************************
    ### Obstacle bounce-back no-slip b.c.
    def bounce_back_obstacle(self, g):

        for k in range(len(self.obstacle)):
            i = self.obstacle[k,0]
            j = self.obstacle[k,1]
            for q in range(1,self.q):
                qb        = self.ns[q]
                dc        = self.c[q, :]
                ic        = self.c[qb,:]
                ii        = i + dc[0]
                jj        = j + dc[1]
                w         = self.lattice[jj,ii]
                if (not w): g[q,j,i] = g[qb,j,i]

        # for q in range(self.q):
        #     g[q,self.lattice] = g[self.ns[q], self.lattice]

    ### ************************************************
    ### Output 2D flow view
    def output_view(self, it, freq, u_in):

        v = self.u.copy()
        v[:,self.obstacle[:,1],self.obstacle[:,0]] = 10.0

        if (it%freq==0):
            plt.clf()
            plt.imshow(np.sqrt(v[0]**2+v[1]**2),
                       cmap = 'RdBu',
                       vmin = 0.0,
                       vmax = u_in,
                       interpolation = 'none')
            filename = self.png_dir+'vel_'+str(self.output_it)+'.png'
            plt.axis('off')
            plt.savefig(filename,
                        dpi=self.dpi)
            self.trim_white(filename)
            self.output_it += 1

    ### ************************************************
    ### Zou-He inlet b.c.
    def zou_he_inlet(self, g, g_eq, rho, u):

        g[1,:,0] = g_eq[1,:,0] + g[2,:,0] - g_eq[2,:,0]
        g[5,:,0] = 0.5*(rho[:,0]*u[0,:,0] - g[1,:,0] + g[2,:,0] \
                 - g[3,:,0] + g[4,:,0] + 2.0*g[6,:,0])
        g[8,:,0] = g[3,:,0] - g[4,:,0] + g[5,:,0] \
                 - g[6,:,0] + g[7,:,0]

    ### ************************************************
    ### Zou-He outlet b.c.
    def zou_he_outlet(self, g, g_eq, rho, u):

        lx             = self.nx-1
        g[2,:,lx] = g_eq[2,:,lx] + g[1,:,lx] \
                       - g_eq[1,:,lx]
        g[6,:,lx] =-0.5*( rho[:,lx]*u[0,:,lx]   \
                             - g[1,:,lx] + g[2,:,lx] \
                             - g[3,:,lx] + g[4,:,lx] \
                             - 2.0*g[5,:,lx])
        g[7,:,lx] =-g[3,:,lx] + g[4,:,lx] \
                       - g[5,:,lx] + g[6,:,lx] \
                       + g[8,:,lx]

    ### ************************************************
    ### Zou-He no-slip top wall b.c.
    def zou_he_top_wall(self, g):

        self.g[4,0,:] = self.g_eq[4,0,:] + self.g[3,0,:] - self.g_eq[3,0,:]
        self.g[6,0,:] = 0.5*(  self.g[1,0,:] + self.g[3,0,:] \
                             - self.g[2,0,:] - self.g[4,0,:] \
                             + 2.0*self.g[5,0,:])
        self.g[8,0,:] = self.g[3,0,:] - self.g[4,0,:] \
                      + self.g[5,0,:] - self.g[6,0,:] \
                      + self.g[7,0,:]

    ### ************************************************
    ### Zou-He no-slip bottom wall b.c.
    def zou_he_bottom_wall(self, g):

        ly             = self.ny-1

        self.g[3,ly,:] = self.g_eq[3,ly,:] + self.g[4,ly,:] \
                       - self.g_eq[4,ly,:]
        self.g[5,ly,:] = 0.5*(2.0*self.g[6,ly,:] - self.g[1,ly,:] \
                                - self.g[3,ly,:] + self.g[2,ly,:] \
                                + self.g[4,ly,:])
        self.g[7,ly,:] = self.g[8,ly,:] - self.g[3,ly,:] \
                       + self.g[4,ly,:] - self.g[5,ly,:] \
                       + self.g[6,ly,:]

    ### ************************************************
    ### Stream distribution
    def stream(self):

        lx                        = self.nx - 1
        ly                        = self.ny - 1
        self.g[0, :, :]           = self.g_up[0, :, :]           # center
        self.g[1, 1:ly-1, 1:lx]   = self.g_up[1, 1:ly-1, 0:lx-1] # +x
        self.g[2, 1:ly-1, 0:lx-1] = self.g_up[2, 1:ly-1, 1:lx]   # -x
        self.g[3, 0:ly-1, 1:lx-1] = self.g_up[3, 1:ly,   1:lx-1] # +y
        self.g[4, 1:ly,   1:lx-1] = self.g_up[4, 0:ly-1, 1:lx-1] # -y
        self.g[5, 0:ly-1, 1:lx]   = self.g_up[5, 1:ly,   0:lx-1] # +x+y
        self.g[6, 1:ly,   0:lx-1] = self.g_up[6, 0:ly-1, 1:lx]   # -x-y
        self.g[7, 0:ly-1, 0:lx-1] = self.g_up[7, 1:ly,   1:lx]   # -x+y
        self.g[8, 1:ly,   1:lx]   = self.g_up[8, 0:ly-1, 0:lx-1] # +x-y

    ### ************************************************
    ### Compute equilibrium state
    def equilibrium(self, g, rho, u):

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
    ### Generate lattice
    def generate(self, polygon):

        # Because we loop on the lattice left-right and top-down,
        # we need to flip the polygon up-down
        poly         = polygon.copy()
        poly[:,1]   *= -1.0
        self.polygon = poly

        # Compute polygon bnds
        poly_bnds    = np.zeros((4))
        poly_bnds[0] = np.amin(poly[:,0])
        poly_bnds[1] = np.amax(poly[:,0])
        poly_bnds[2] = np.amin(poly[:,1])
        poly_bnds[3] = np.amax(poly[:,1])

        # Declare lattice arrays
        self.lattice  = np.zeros((self.ny, self.nx), dtype=bool)
        obstacle      = np.empty((0,2),              dtype=int)
        self.obstacle = np.empty((0,2),              dtype=int)

        # Fill lattice
        bar = progress.bar.Bar('Generating...', max=self.nx*self.ny)
        for i in range(self.nx):
            for j in range(self.ny):
                pt     = self.lattice_coords(j, i)
                inside = False

                # Check if pt is inside polygon bbox
                if ((pt[0] > poly_bnds[0]) and (pt[0] < poly_bnds[1])):
                    if ((pt[1] > poly_bnds[2]) and (pt[1] < poly_bnds[3])):
                        inside = self.is_inside(poly, pt)

                        if (inside):
                            obstacle = np.append(obstacle,
                                                 np.array([[i,j]]),
                                                 axis=0)

                # Fill lattice
                self.lattice[j,i] = inside

                bar.next()
        bar.finish()

        print('Found '+str(obstacle.shape[0])+' locations in obstacle')

        # Re-process obstacle to keep boundary
        for k in range(obstacle.shape[0]):
            i   = obstacle[k,0]
            j   = obstacle[k,1]
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if (not self.lattice[j+dj,i+di]):
                        self.obstacle = np.append(self.obstacle,
                                                  np.array([[i,j]]),
                                                  axis=0)

        # Printings
        print('Found '+str(self.obstacle.shape[0])+' on obstacle boundary')

    ### ************************************************
    ### Get lattice coordinates from integers
    def lattice_coords(self, j, i):

        # Compute and return the coordinates of the lattice node (i,j)
        dx = (self.xmax - self.xmin)/(self.nx - 1)
        dy = (self.ymax - self.ymin)/(self.ny - 1)
        x  = self.xmin + i*dx
        y  = self.ymin + (self.ny-j)*dy

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
        filename = self.output_dir+self.name+'.png'

        plt.axis('off')
        plt.imshow(self.lattice,
                   cmap=mplt.cm.inferno,
                   interpolation='none')
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        plt.close()
        self.trim_white(filename)

    ### ************************************************
    ### Crop white background from image
    def trim_white(self, filename):

        # Trim using PIL
        im   = Image.open(filename)
        bg   = Image.new(im.mode, im.size, (255,255,255))
        diff = PIL.ImageChops.difference(im, bg)
        bbox = diff.getbbox()
        cp   = im.crop(bbox)
        cp.save(filename)

    ### ************************************************
    ### Poiseuille flow
    def poiseuille(self, d, y, x):

        H = self.ny
        return (1.0-d)*4.0*y*(H-y)/H**2
