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
                 tau_p_lbm, tau_m_lbm,
                 Cx,         Ct, Cs,
                 Cr,         Cu,
                 Cf,         dx, dt,
                 output_dir, dpi):

        self.name       = name
        self.xmin       = xmin
        self.xmax       = xmax
        self.ymin       = ymin
        self.ymax       = ymax
        self.nx         = nx
        self.ny         = ny
        self.q          = 9
        self.tau_p_lbm = tau_p_lbm
        self.tau_m_lbm = tau_m_lbm
        self.dx         = dx
        self.dt         = dt
        self.Cx         = Cx
        self.Ct         = Ct
        self.Cs         = Cs
        self.Cr         = Cr
        self.Cu         = Cu
        self.Cf         = Cf
        self.output_dir = output_dir
        self.png_dir    = self.output_dir+'./png/'
        self.output_it  = 0
        self.dpi        = dpi
        self.lx         = self.nx - 1
        self.ly         = self.ny - 1

        if (not os.path.exists(self.output_dir)):
            os.makedirs(self.output_dir)
        if (not os.path.exists(self.png_dir)):
            os.makedirs(self.png_dir)

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
        # Lattice array is oriented as follows :
        # +x     = left-right
        # +y     = bottom-top
        # origin = bottom left
        self.g       = np.zeros((self.q,  self.nx, self.ny))
        #self.g_s     = np.zeros((self.q,  self.nx, self.ny))
        self.g_eq    = np.zeros((self.q,  self.nx, self.ny))
        #self.g_m     = np.zeros((self.q,  self.nx, self.ny))
        #self.g_p     = np.zeros((self.q,  self.nx, self.ny))
        self.g_up    = np.zeros((self.q,  self.nx, self.ny))

        # Lattice array is oriented as follows :
        # +x     = left-right
        # +y     = bottom-top
        # origin = bottom left
        self.lattice = np.zeros((self.nx, self.ny), dtype=bool)

    ### ************************************************
    ### Solve LBM
    def solve(self, it_max, u_lbm, rho_lbm,
                    R_ref,  U_ref, L_ref,
                    freq):

        self.rho     = np.ones ((   self.nx, self.ny))
        self.u       = np.zeros((2, self.nx, self.ny))

        # Input profiles
        self.input_velocity(u_lbm)
        self.u[:,0,:] = self.u_in[:,:]
        self.rho     *= rho_lbm
        #self.u[:,self.lattice] = 0.0

        self.equilibrium()
        self.g = self.g_eq

        # Solve
        bar = progress.bar.Bar('Solving...', max=it_max)
        for it in range(it_max+1):

            # Drag and lift
            #self.drag_lift(it, R_ref, U_ref, L_ref)

            # Compute macroscopic fields
            self.macro()

            # Compute equilibrium state
            self.equilibrium()

            # Output view
            self.output_view(it, freq, u_lbm)

            # Collisions
            self.trt_collisions()

            #self.g_s = self.g.copy()
            #self.bounce_back_obstacle(self.g, self.g_up)

            # Streaming
            self.stream()

            # Boundary conditions
            self.zou_he_top_wall()
            self.zou_he_bottom_wall()
            self.zou_he_inlet()
            self.zou_he_outlet(rho_lbm)
            self.zou_he_bottom_left_corner()
            self.zou_he_top_left_corner()
            self.zou_he_top_right_corner()
            self.zou_he_bottom_right_corner()

            # Increment bar
            bar.next()

        # End bar
        bar.finish()

    ### ************************************************
    ### TRT collision operator
    def trt_collisions(self):

        #self.g_up = self.g - (1.0/self.tau_p_lbm)*(self.g - self.g_eq)

        # Compute g_p = g_p - g_eq_p
        #     and g_m = g_m - g_eq_m
        self.g_p = 0.5*(self.g[:,:,:]    + self.g[self.ns[:],:,:] \
                     - (self.g_eq[:,:,:] + self.g_eq[self.ns[:],:,:]))
        self.g_m = 0.5*(self.g[:,:,:]    - self.g[self.ns[:],:,:] \
                     - (self.g_eq[:,:,:] - self.g_eq[self.ns[:],:,:]))
        self.g_m[0,:,:] += 0.5*(self.g[0,:,:] - self.g_eq[0,:,:])

        # Compute collisions
        self.g_up = self.g - (1.0/self.tau_p_lbm)*self.g_p \
                          - (1.0/self.tau_m_lbm)*self.g_m

    ### ************************************************
    ### Compute input velocity
    def input_velocity(self, u_lbm):

        self.u_in = np.zeros((2, self.ny))

        for j in range(self.ny):
            pt             = self.lattice_coords(0, j)
            self.u_in[:,j] = u_lbm*self.poiseuille(pt)

    ### ************************************************
    ### Compute drag and lift
    def drag_lift(self, it, R_ref, U_ref, L_ref):

        # Initialize
        force = np.zeros((2))

        # Loop over obstacle array
        for k in range(len(self.boundary)):
            i = self.boundary[k,0]
            j = self.boundary[k,1]

            for q in range(1,self.q):
                qb        = self.ns[q]
                dc        = self.c[q, :]
                ic        = self.c[qb,:]
                ii        = i + dc[0]
                jj        = j + dc[1]
                w         = self.lattice[ii,jj]
                df        =-self.g[qb,i,j]*ic[:] + self.g_up[q,i,j]*dc[:]
                force[:] += w*df

        # Normalize coefficient
        force *= self.Cf
        force *= 1.0/(0.5*R_ref*L_ref*U_ref**2)

        # Write to file
        filename = self.output_dir+'drag_lift'
        with open(filename, 'a') as f:
            f.write('{} {} {}\n'.format(it*self.dt, force[0], force[1]))

    ### ************************************************
    ### Obstacle halfway bounce-back no-slip b.c.
    def bounce_back_obstacle(self, g, g_up):

        for k in range(len(self.boundary)):
            i = self.boundary[k,0]
            j = self.boundary[k,1]
            for q in range(1,self.q):
                qb = self.ns[q]
                dc = self.c[q, :]
                ic = self.c[qb,:]
                ii = i + dc[0]
                jj = j + dc[1]
                w  = self.lattice[ii,jj]

                # Apply if neighbor is solid
                if w: self.g[qb,i,j] = self.g[q,i,j]

    ### ************************************************
    ### Zou-He inlet b.c.
    def zou_he_inlet(self):

        self.u[0,0,:] = self.u_in[0,:]
        self.u[1,0,:] = self.u_in[1,:]

        self.rho[0,:] = (self.g[0,0,:] +
                         self.g[3,0,:] +
                         self.g[4,0,:] +
                     2.0*self.g[2,0,:] +
                     2.0*self.g[6,0,:] +
                     2.0*self.g[7,0,:] )/(1.0 - self.u_in[0,:])

        self.g[1,0,:] = (self.g_eq[1,0,:] +
                         self.g[2,0,:]    -
                         self.g_eq[2,0,:] )

        self.g[5,0,:] =-0.5*(self.g[1,0,:] -
                             self.g[2,0,:] +
                             self.g[3,0,:] -
                             self.g[4,0,:] -
                         2.0*self.g[6,0,:] -
                             self.rho[0,:]*self.u_in[0,:] -
                             self.rho[0,:]*self.u_in[1,:])

        self.g[8,0,:] =-0.5*(self.g[1,0,:] -
                             self.g[2,0,:] -
                             self.g[3,0,:] +
                             self.g[4,0,:] -
                         2.0*self.g[7,0,:] -
                             self.rho[0,:]*self.u_in[0,:] +
                             self.rho[0,:]*self.u_in[1,:])

    ### ************************************************
    ### Zou-He outlet b.c.
    def zou_he_outlet(self, rho_lbm):

        lx = self.lx

        self.rho[lx,:] = rho_lbm
        self.u[1,lx,:] = 0.0

        self.u[0,lx,:] = (self.g[0,lx,:] +
                          self.g[3,lx,:] +
                          self.g[4,lx,:] +
                      2.0*self.g[1,lx,:] +
                      2.0*self.g[5,lx,:] +
                      2.0*self.g[8,lx,:])/self.rho[lx,:] - 1.0

        self.g[2,lx,:] = (self.g_eq[2,lx,:] +
                          self.g[1,lx,:]    -
                          self.g_eq[1,lx,:])

        self.g[6,lx,:] = 0.5*(self.g[1,lx,:] -
                              self.g[2,lx,:] +
                              self.g[3,lx,:] -
                              self.g[4,lx,:] +
                          2.0*self.g[5,lx,:] -
                              self.rho[lx,:]*self.u[0,lx,:] -
                              self.rho[lx,:]*self.u[1,lx,:])

        self.g[7,lx,:] = 0.5*(self.g[1,lx,:] -
                              self.g[2,lx,:] -
                              self.g[3,lx,:] +
                              self.g[4,lx,:] +
                          2.0*self.g[8,lx,:] -
                              self.rho[lx,:]*self.u[0,lx,:] +
                              self.rho[lx,:]*self.u[1,lx,:])

    ### ************************************************
    ### Zou-He no-slip top wall b.c.
    def zou_he_top_wall(self):

        ly = self.ny - 1

        u_x = 0.0
        u_y = 0.0
        self.u[0,:,ly] = u_x
        self.u[1,:,ly] = u_y

        self.rho[:,ly] = (self.g[0,:,ly] +
                          self.g[1,:,ly] +
                          self.g[2,:,ly] +
                      2.0*self.g[3,:,ly] +
                      2.0*self.g[5,:,ly] +
                      2.0*self.g[7,:,ly] )/(1.0 + u_y)

        self.g[4,:,ly] = (self.g_eq[4,:,ly] +
                          self.g[3,:,ly]    -
                          self.g_eq[3,:,ly] )

        self.g[6,:,ly] = 0.5*(self.g[1,:,ly] -
                              self.g[2,:,ly] +
                              self.g[3,:,ly] -
                              self.g[4,:,ly] +
                          2.0*self.g[5,:,ly] -
                              self.rho[:,ly]*u_x -
                              self.rho[:,ly]*u_y)

        self.g[8,:,ly] =-0.5*(self.g[1,:,ly] -
                              self.g[2,:,ly] -
                              self.g[3,:,ly] +
                              self.g[4,:,ly] -
                          2.0*self.g[7,:,ly] -
                              self.rho[:,ly]*u_x +
                              self.rho[:,ly]*u_y)


    ### ************************************************
    ### Zou-He no-slip bottom wall b.c.
    def zou_he_bottom_wall(self):

        u_x = 0.0
        u_y = 0.0
        self.u[0,:,0] = u_x
        self.u[1,:,0] = u_y

        self.rho[:,0] = (self.g[0,:,0] +
                         self.g[1,:,0] +
                         self.g[2,:,0] +
                     2.0*self.g[4,:,0] +
                     2.0*self.g[6,:,0] +
                     2.0*self.g[8,:,0] )/(1.0 - u_x)

        self.g[3,:,0] = (self.g_eq[3,:,0] +
                         self.g[4,:,0]    -
                         self.g_eq[4,:,0] )

        self.g[5,:,0] =-0.5*(self.g[1,:,0] -
                             self.g[2,:,0] +
                             self.g[3,:,0] -
                             self.g[4,:,0] -
                         2.0*self.g[6,:,0] -
                             self.rho[:,0]*u_x -
                             self.rho[:,0]*u_y)

        self.g[7,:,0] = 0.5*(self.g[1,:,0] -
                             self.g[2,:,0] -
                             self.g[3,:,0] +
                             self.g[4,:,0] +
                         2.0*self.g[8,:,0] -
                             self.rho[:,0]*u_x +
                             self.rho[:,0]*u_y)

    ### ************************************************
    ### Zou-He bottom left corner
    def zou_he_bottom_left_corner(self):

        u_x = 0.0
        u_y = 0.0
        self.u[0,0,0] = u_x
        self.u[1,0,0] = u_y

        self.rho[0,0] = self.rho[1,0]

        self.g[1,0,0] = (self.g_eq[1,0,0] +
                         self.g[2,0,0]    -
                         self.g_eq[2,0,0] )

        self.g[3,0,0] = (self.g_eq[3,0,0] +
                         self.g[4,0,0]    -
                         self.g_eq[4,0,0] )

        self.g[5,0,0] = (self.g_eq[5,0,0] +
                         self.g[6,0,0]    -
                         self.g_eq[6,0,0] )

        self.g[7,0,0] = 0.0
        self.g[8,0,0] = 0.0

        self.g[0,0,0] = (self.rho[0,0] -
                         self.g[1,0,0] -
                         self.g[2,0,0] -
                         self.g[3,0,0] -
                         self.g[4,0,0] -
                         self.g[5,0,0] -
                         self.g[6,0,0] -
                         self.g[7,0,0] -
                         self.g[8,0,0] )

    ### ************************************************
    ### Zou-He top left corner
    def zou_he_top_left_corner(self):

        ly = self.ly

        u_x = 0.0
        u_y = 0.0
        self.u[0,0,ly] = u_x
        self.u[1,0,ly] = u_y

        self.rho[0,ly] = self.rho[1,ly]

        self.g[1,0,ly] = (self.g_eq[1,0,ly] +
                          self.g[2,0,ly]    -
                          self.g_eq[2,0,ly] )

        self.g[4,0,ly] = (self.g_eq[4,0,ly] +
                          self.g[3,0,ly]    -
                          self.g_eq[3,0,ly] )

        self.g[8,0,ly] = (self.g_eq[8,0,ly] +
                          self.g[7,0,ly]    -
                          self.g_eq[7,0,ly] )

        self.g[5,0,ly] = 0.0
        self.g[6,0,ly] = 0.0

        self.g[0,0,ly] = (self.rho[0,ly] -
                          self.g[1,0,ly] -
                          self.g[2,0,ly] -
                          self.g[3,0,ly] -
                          self.g[4,0,ly] -
                          self.g[5,0,ly] -
                          self.g[6,0,ly] -
                          self.g[7,0,ly] -
                          self.g[8,0,ly] )

    ### ************************************************
    ### Zou-He top right corner
    def zou_he_top_right_corner(self):

        lx = self.lx
        ly = self.ly

        #u_y = 0.0
        self.u[0,lx,ly] = self.u[0,lx-1,ly]
        self.u[1,lx,ly] = self.u[1,lx-1,ly]
        self.rho[lx,ly] = self.rho[lx-1,ly]

        self.g[2,lx,ly] = (self.g_eq[2,lx,ly] +
                          self.g[1,lx,ly]    -
                          self.g_eq[1,lx,ly] )

        self.g[4,lx,ly] = (self.g_eq[4,lx,ly] +
                          self.g[3,lx,ly]    -
                          self.g_eq[3,lx,ly] )

        self.g[6,lx,ly] = (self.g_eq[6,lx,ly] +
                          self.g[5,lx,ly]    -
                          self.g_eq[5,lx,ly] )

        self.g[7,lx,ly] = 0.0
        self.g[8,lx,ly] = 0.0

        self.g[0,lx,ly] = (self.rho[lx,ly] -
                           self.g[1,lx,ly] -
                           self.g[2,lx,ly] -
                           self.g[3,lx,ly] -
                           self.g[4,lx,ly] -
                           self.g[5,lx,ly] -
                           self.g[6,lx,ly] -
                           self.g[7,lx,ly] -
                           self.g[8,lx,ly] )

    ### ************************************************
    ### Zou-He bottom right corner
    def zou_he_bottom_right_corner(self):

        lx = self.lx

        #u_y = 0.0
        self.u[0,lx,0] = self.u[0,lx-1,0]
        self.u[1,lx,0] = self.u[1,lx-1,0]
        self.rho[lx,0] = self.rho[lx-1,0]

        self.g[2,lx,0] = (self.g_eq[2,lx,0] +
                          self.g[1,lx,0]    -
                          self.g_eq[1,lx,0] )

        self.g[3,lx,0] = (self.g_eq[3,lx,0] +
                          self.g[4,lx,0]    -
                          self.g_eq[4,lx,0] )

        self.g[7,lx,0] = (self.g_eq[7,lx,0] +
                          self.g[8,lx,0]    -
                          self.g_eq[8,lx,0] )

        self.g[5,lx,0] = 0.0
        self.g[6,lx,0] = 0.0

        self.g[0,lx,0] = (self.rho[lx,0] -
                          self.g[1,lx,0] -
                          self.g[2,lx,0] -
                          self.g[3,lx,0] -
                          self.g[4,lx,0] -
                          self.g[5,lx,0] -
                          self.g[6,lx,0] -
                          self.g[7,lx,0] -
                          self.g[8,lx,0] )

    ### ************************************************
    ### Stream distribution
    def stream(self):

        for q in range(self.q):
            self.g[q,:,:] = np.roll(
                            np.roll(
                                self.g_up[q,:,:],self.c[q,0],axis=0),
                                                 self.c[q,1],axis=1)


        # lx                        = self.lx
        # ly                        = self.ly
        # self.g[0, :, :]           = self.g_up[0, :, :]           # center
        # self.g[1, 1:lx,   0:ly]   = self.g_up[1, 0:lx-1, 0:ly]   # +x
        # self.g[2, 0:lx-1, 0:ly]   = self.g_up[2, 1:lx,   0:ly]   # -x
        # self.g[3, 0:lx,   1:ly]   = self.g_up[3, 0:lx,   0:ly-1] # +y
        # self.g[4, 0:lx,   0:ly-1] = self.g_up[4, 0:lx,   1:ly]   # -y
        # self.g[5, 1:lx,   1:ly]   = self.g_up[5, 0:lx-1, 0:ly-1] # +x+y
        # self.g[6, 0:lx-1, 0:ly-1] = self.g_up[6, 1:lx,   1:ly]   # -x-y
        # self.g[7, 0:lx-1, 1:ly]   = self.g_up[7, 1:lx,   0:ly-1] # -x+y
        # self.g[8, 1:lx,   0:ly-1] = self.g_up[8, 0:lx-1, 1:ly]   # +x-y

        # self.g[:,self.lattice]    = self.g_s[:,self.lattice]

    ### ************************************************
    ### Compute equilibrium state
    def equilibrium(self):

        # Compute velocity term
        v = (3.0/2.0)*(self.u[0,:,:]**2 + self.u[1,:,:]**2)

        # Compute equilibrium
        for q in range(self.q):
            t                 = 3.0*(self.u[0,:,:]*self.c[q,0] +
                                     self.u[1,:,:]*self.c[q,1])
            self.g_eq[q,:,:]  = (1.0 + t + 0.5*t**2 - v)
            self.g_eq[q,:,:] *= self.rho[:,:]*self.w[q]

    ### ************************************************
    ### Compute macroscopic fields
    def macro(self):

        self.macro_density()
        self.macro_velocity()

    ### ************************************************
    ### Compute macroscopic density
    def macro_density(self):

        # Compute density
        self.rho = np.sum(self.g, axis=0)

    ### ************************************************
    ### Compute macroscopic velocity
    def macro_velocity(self):

        # Compute velocity
        self.u[:,:,:] = 0.0

        for q in range(self.q):
            self.u[0,:,:] += self.c[q,0]*self.g[q,:,:]
            self.u[1,:,:] += self.c[q,1]*self.g[q,:,:]

        self.u[0,:,:] /= self.rho[:,:]
        self.u[1,:,:] /= self.rho[:,:]

    ### ************************************************
    ### Output 2D flow view
    def output_view(self, it, freq, u_in):

        v = self.u.copy()
        #v[:,self.boundary[:,1],self.boundary[:,0]] = 10.0

        if (it%freq==0):
            plt.clf()
            plt.imshow(np.rot90(v[0]),
                       cmap = 'coolwarm',
                       vmin = 0.0,
                       vmax = u_in,
                       interpolation = 'nearest')
            filename = self.png_dir+'vel_'+str(self.output_it)+'.png'
            plt.axis('off')
            plt.savefig(filename,
                        dpi=self.dpi)
            self.trim_white(filename)
            self.output_it += 1

    ### ************************************************
    ### Add obstacle
    def add_obstacle(self, polygon):

        # Because we loop on the lattice left-right and top-down,
        # we need to flip the polygon up-down
        poly         = polygon.copy()
        #poly[:,1]   *= -1.0
        self.polygon = poly

        # Compute polygon bnds
        poly_bnds    = np.zeros((4))
        poly_bnds[0] = np.amin(poly[:,0])
        poly_bnds[1] = np.amax(poly[:,0])
        poly_bnds[2] = np.amin(poly[:,1])
        poly_bnds[3] = np.amax(poly[:,1])

        # Declare lattice arrays
        obstacle      = np.empty((0,2), dtype=int)
        self.boundary = np.empty((0,2), dtype=int)

        # Fill lattice
        bar = progress.bar.Bar('Generating...', max=self.nx*self.ny)
        for i in range(self.nx):
            for j in range(self.ny):
                pt = self.lattice_coords(i, j)

                # Check if pt is inside polygon bbox
                if ((pt[0] > poly_bnds[0]) and
                    (pt[0] < poly_bnds[1]) and
                    (pt[1] > poly_bnds[2]) and
                    (pt[1] < poly_bnds[3])):

                    if (self.is_inside(poly, pt)):
                        obstacle = np.append(obstacle,
                                             np.array([[i,j]]), axis=0)
                        self.lattice[i,j] = True

                bar.next()
        bar.finish()

        print('Found '+str(obstacle.shape[0])+' locations in obstacle')

        # Check area of obstacle
        self.area = 0.0

        for i in range(self.nx):
            for j in range(self.ny):
                if (self.lattice[i,j]): self.area += self.dx**2

        print('Obstacle area: '+str(self.area))

        # Re-process obstacle to keep first solid nodes only
        for k in range(obstacle.shape[0]):
            i = obstacle[k,0]
            j = obstacle[k,1]
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if (not self.lattice[i+di,j+dj]):
                        self.boundary = np.append(self.boundary,
                                                  np.array([[i+di,j+dj]]),
                                                  axis=0)

        # Printings
        print('Found '+str(self.boundary.shape[0])+' on obstacle boundary')

    ### ************************************************
    ### Get lattice coordinates from integers
    def lattice_coords(self, i, j):

        # Compute and return the coordinates of the lattice node (i,j)
        dx = (self.xmax - self.xmin)/(self.nx - 1)
        dy = (self.ymax - self.ymin)/(self.ny - 1)
        x  = i*dx
        y  = j*dy

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

        # Add obstacle border
        lat = self.lattice.copy()
        lat = lat.astype(float)
        #lat[self.boundary[:,1],self.boundary[:,0]] += 2.0

        # Plot and save image of lattice
        filename = self.output_dir+self.name+'.png'

        plt.axis('off')
        plt.imshow(np.transpose(lat),
                   cmap = mplt.cm.inferno,
                   vmin = 0.0,
                   vmax = 2.0,
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
    def poiseuille(self, pt):

        x    = pt[0]
        y    = pt[1]
        H    = self.ymax - self.ymin
        u    = np.zeros(2)
        u[0] = 4.0*y*(H-y)/H**2

        return u
