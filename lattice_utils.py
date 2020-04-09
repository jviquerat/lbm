# Generic imports
import os
import PIL
import math
import random
import numpy             as np
import matplotlib        as mplt
import matplotlib.pyplot as plt
from   PIL               import Image
from   datetime          import datetime
from   numba             import jit

# Custom imports
from   buff              import *

### ************************************************
### Class defining an obstacle in the lattice
class Obstacle:
    ### ************************************************
    ### Constructor
    def __init__(self, polygon, area, boundary, ibb):

        self.polygon  = polygon
        self.area     = area
        self.boundary = boundary
        self.ibb      = ibb

### ************************************************
### Class defining lattice object
class Lattice:
    ### ************************************************
    ### Constructor
    def __init__(self, *args, **kwargs):

        # Input parameters
        self.name       = kwargs.get('name',     'lattice'                  )
        self.x_min      = kwargs.get('x_min',     0.0                       )
        self.x_max      = kwargs.get('x_max',     1.0                       )
        self.y_min      = kwargs.get('y_min',     0.0                       )
        self.y_max      = kwargs.get('y_max',     1.0                       )
        self.nx         = kwargs.get('nx',        100                       )
        self.ny         = kwargs.get('ny',        self.nx                   )
        self.tau_lbm    = kwargs.get('tau_lbm',   1.0                       )
        self.dx         = kwargs.get('dx',        1.0                       )
        self.dt         = kwargs.get('dt',        1.0                       )
        self.Cx         = kwargs.get('Cx',        self.dx                   )
        self.Ct         = kwargs.get('Ct',        self.dt                   )
        self.Cr         = kwargs.get('Cr',        1.0                       )
        self.Cn         = kwargs.get('Cn',        self.Cx**2/self.Ct        )
        self.Cu         = kwargs.get('Cu',        self.Cx/self.Ct           )
        self.Cf         = kwargs.get('Cf',        self.Cr*self.Cx**2/self.Ct)
        self.dpi        = kwargs.get('dpi',       100                       )
        self.u_lbm      = kwargs.get('u_lbm',     0.05                      )
        self.L_lbm      = kwargs.get('L_lbm',     100.0                     )
        self.nu_lbm     = kwargs.get('nu_lbm',    0.01                      )
        self.Re_lbm     = kwargs.get('Re_lbm',    100.0                     )
        self.rho_lbm    = kwargs.get('rho_lbm',   1.0                       )
        self.IBB        = kwargs.get('IBB',       False                     )
        self.stop       = kwargs.get('stop',      'it'                      )
        self.t_max      = kwargs.get('t_max',     1.0                       )
        self.it_max     = kwargs.get('it_max',    1000                      )
        self.obs_cv_ct  = kwargs.get('obs_cv_ct', 1.0e-2                    )
        self.obs_cv_nb  = kwargs.get('obs_cv_nb', 1000                      )

        # Other parameters
        self.output_it  = 0
        self.lx         = self.nx - 1
        self.ly         = self.ny - 1
        self.q          = 9
        self.Cs         = 1.0/math.sqrt(3.0)

        # Output dirs
        time             = datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
        self.results_dir = './results/'
        self.output_dir  = self.results_dir+str(time)+'/'
        self.png_dir     = self.output_dir+'./png/'

        if (not os.path.exists(self.results_dir)):
            os.makedirs(self.results_dir)
        if (not os.path.exists(self.output_dir)):
            os.makedirs(self.output_dir)
        if (not os.path.exists(self.png_dir)):
            os.makedirs(self.png_dir)

        # TRT parameters
        self.tau_p_lbm  = self.tau_lbm
        self.lambda_trt = 1.0/4.0 # Best for stability
        self.tau_m_lbm  = self.lambda_trt/(self.tau_p_lbm - 0.5) + 0.5
        self.om_p_lbm   = 1.0/self.tau_p_lbm
        self.om_m_lbm   = 1.0/self.tau_m_lbm
        self.om_lbm     = 1.0/self.tau_lbm

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

        # Array for bounce-back
        self.ns = np.array([0,2,1,4,3,6,5,8,7])

        # Density arrays
        self.g       = np.zeros((self.q,  self.nx, self.ny))
        self.g_eq    = np.zeros((self.q,  self.nx, self.ny))
        self.g_up    = np.zeros((self.q,  self.nx, self.ny))

        # Lattice array is oriented as follows :
        # +x     = left-right
        # +y     = bottom-top
        # origin = bottom left
        self.lattice = np.zeros((self.nx, self.ny), dtype=bool)

        # Physical fields
        self.rho     = np.ones ((   self.nx, self.ny))
        self.u       = np.zeros((2, self.nx, self.ny))

        # Obstacles
        self.obstacles = []

        # Iterating and stopping
        self.it        = 0
        self.compute   = True
        self.drag_buff = Buff('drag',
                              self.dt,
                              self.obs_cv_ct,
                              self.obs_cv_nb,
                              self.output_dir)
        self.lift_buff = Buff('lift',
                              self.dt,
                              self.obs_cv_ct,
                              self.obs_cv_nb,
                              self.output_dir)

        # Printings
        print('### LBM solver ###')
        print('# u_lbm      = '+str(self.u_lbm))
        print('# L_lbm      = '+str(self.L_lbm))
        print('# nu_lbm     = '+str(self.nu_lbm))
        print('# Re_lbm     = '+str(self.Re_lbm))
        print('# tau_p_lbm  = '+str(self.tau_p_lbm))
        print('# tau_m_lbm  = '+str(self.tau_m_lbm))
        print('# dt         = '+str(self.dt))
        print('# dx         = '+str(self.dx))
        print('# nx         = '+str(self.nx))
        print('# ny         = '+str(self.ny))

    ### ************************************************
    ### Compute macroscopic fields
    def macro(self):

        # Compute density
        self.rho[:,:] = np.sum(self.g[:,:,:], axis=0)

        # Compute velocity
        self.u[0,:,:] = np.tensordot(self.c[:,0],
                                     self.g[:,:,:],
                                     axes=(0,0))
        self.u[1,:,:] = np.tensordot(self.c[:,1],
                                     self.g[:,:,:],
                                     axes=(0,0))

        self.u[0,:,:] /= self.rho[:,:]
        self.u[1,:,:] /= self.rho[:,:]

    ### ************************************************
    ### Compute equilibrium state
    def equilibrium(self):

        nb_equilibrium(self.u, self.c, self.w, self.rho, self.g_eq)

        # # Compute velocity term
        # v = (3.0/2.0)*(self.u[0,:,:]**2 + self.u[1,:,:]**2)

        # # Compute equilibrium
        # for q in range(self.q):
        #     t                 = 3.0*(self.u[0,:,:]*self.c[q,0] +
        #                              self.u[1,:,:]*self.c[q,1])
        #     self.g_eq[q,:,:]  = (1.0 + t + 0.5*t**2 - v)
        #     self.g_eq[q,:,:] *= self.rho[:,:]*self.w[q]

    ### ************************************************
    ### Collision and streaming
    def collision_stream(self):

        nb_col_str(self.g, self.g_eq, self.g_up,
                   self.om_p_lbm, self.om_m_lbm,
                   self.c, self.ns,
                   self.nx, self.ny,
                   self.lx, self.ly)

        # # Take care of q=0 first
        # self.g_up[0,:,:] = (self.g[0,:,:]
        #                  -  self.om_p_lbm*(self.g[0,:,:] - self.g_eq[0,:,:]))
        # self.g   [0,:,:] =  self.g_up[0,:,:]

        # # Collide other indices
        # for q in range(1,self.q):
        #     qb = self.ns[q]

        #     self.g_up[q,:,:] = (self.g   [q,:,:]
        #    - self.om_p_lbm*0.5*(self.g   [q,:,:]
        #                       + self.g   [qb,:,:]
        #                       - self.g_eq[q,:,:]
        #                       - self.g_eq[qb,:,:])
        #    - self.om_m_lbm*0.5*(self.g   [q,:,:]
        #                       - self.g   [qb,:,:]
        #                       - self.g_eq[q,:,:]
        #                       + self.g_eq[qb,:,:]))

        # # Stream
        # nx = self.nx
        # ny = self.ny
        # lx = self.lx
        # ly = self.ly

        # self.g[1,1:nx, :  ] = self.g_up[1,0:lx, :  ]
        # self.g[2,0:lx, :  ] = self.g_up[2,1:nx, :  ]
        # self.g[3, :,  1:ny] = self.g_up[3, :,  0:ly]
        # self.g[4, :,  0:ly] = self.g_up[4, :,  1:ny]
        # self.g[5,1:nx,1:ny] = self.g_up[5,0:lx,0:ly]
        # self.g[6,0:lx,0:ly] = self.g_up[6,1:nx,1:ny]
        # self.g[7,0:lx,1:ny] = self.g_up[7,1:nx,0:ly]
        # self.g[8,1:nx,0:ly] = self.g_up[8,0:lx,1:ny]

    ### ************************************************
    ### Compute drag and lift
    def drag_lift(self, obs, R_ref, U_ref, L_ref):

        # Initialize
        fx     = 0.0
        fy     = 0.0

        # Loop over obstacle array
        for k in range(len(self.obstacles[obs].boundary)):
            i    = self.obstacles[obs].boundary[k,0]
            j    = self.obstacles[obs].boundary[k,1]
            q    = self.obstacles[obs].boundary[k,2]
            qb   = self.ns[q]
            cx   = self.c[q,0]
            cy   = self.c[q,1]
            g_up = self.g_up[q,i,j] + self.g[qb,i,j]

            fx += g_up*cx
            fy += g_up*cy

        # Normalize coefficient
        Cx =-2.0*fx/(R_ref*L_ref*U_ref**2)
        Cy =-2.0*fy/(R_ref*L_ref*U_ref**2)

        return Cx, Cy

    ### ************************************************
    ### Handle drag/lift buffers
    def add_buff(self, Cx, Cy, it):

        # Add to buffer and check for convergence
        self.drag_buff.add(Cx)
        self.lift_buff.add(Cy)

        avg_Cx, dcx = self.drag_buff.mv_avg()
        avg_Cy, dcy = self.lift_buff.mv_avg()

        # Write to file
        filename = self.output_dir+'drag_lift'
        with open(filename, 'a') as f:
            f.write('{} {} {} {} {} {} {}\n'.format(it*self.dt,
                                              Cx,     Cy,
                                              avg_Cx, avg_Cy,
                                              dcx,    dcy))

        return avg_Cx, avg_Cy

    ### ************************************************
    ### Obstacle halfway bounce-back no-slip b.c.
    def bounce_back_obstacle(self, obs):

        # Interpolated BB
        if (self.IBB):
            for k in range(len(self.obstacles[obs].boundary)):
                i  = self.obstacles[obs].boundary[k,0]
                j  = self.obstacles[obs].boundary[k,1]
                q  = self.obstacles[obs].boundary[k,2]
                qb = self.ns[q]
                c  = self.c[q,:]
                cb = self.c[qb,:]
                im = i + cb[0]
                jm = j + cb[1]
                imm = i + 2*cb[0]
                jmm = j + 2*cb[1]

                p  = self.obstacles[obs].ibb[k]
                pp = 2.0*p
                if (p < 0.5):
                    self.g[qb,i,j] = (p*(pp+1.0)*self.g_up[q,i,j]
                                   + (1.0+pp)*(1.0-pp)*self.g_up[q,im,jm]
                                   - p*(1.0-pp)*self.g_up[q,imm,jmm])
                else:
                    self.g[qb,i,j] = ((1.0/(p*(pp+1.0)))*self.g_up[q,i,j] +
                                     ((pp-1.0)/p)*self.g_up[qb,i,j] +
                                     ((1.0-pp)/(1.0+pp))*self.g_up[qb,im,jm])

        # Regular BB
        if (not self.IBB):
            for k in range(len(self.obstacles[obs].boundary)):
                i  = self.obstacles[obs].boundary[k,0]
                j  = self.obstacles[obs].boundary[k,1]
                q  = self.obstacles[obs].boundary[k,2]
                qb = self.ns[q]
                c  = self.c[q,:]
                ii = i + c[0]
                jj = j + c[1]

                self.g[qb,i,j] = self.g_up[q,i,j]

        # Set velocity of obstacles
        self.u[:,self.lattice] = 0.0

    ### ************************************************
    ### Zou-He left wall velocity b.c.
    def zou_he_left_wall_velocity(self):

        lx = self.lx
        ly = self.ly

        self.u[0,0,:] = self.u_left[0,:]
        self.u[1,0,:] = self.u_left[1,:]

        self.rho[0,:] = (self.g[0,0,:] +
                         self.g[3,0,:] +
                         self.g[4,0,:] +
                     2.0*self.g[2,0,:] +
                     2.0*self.g[6,0,:] +
                     2.0*self.g[7,0,:] )/(1.0 - self.u[0,0,:])

        self.g[1,0,:] = (self.g[2,0,:] +
                         (2.0/3.0)*self.rho[0,:]*self.u[0,0,:])

        self.g[5,0,:] = (     self.g[6,0,:]  -
                         0.5*(self.g[3,0,:]  -
                              self.g[4,0,:]) +
                    (1.0/6.0)*self.rho[0,:]*self.u[0,0,:] +
                    (1.0/2.0)*self.rho[0,:]*self.u[1,0,:] )

        self.g[8,0,:] = (     self.g[7,0,:]  +
                         0.5*(self.g[3,0,:]  -
                              self.g[4,0,:]) +
                    (1.0/6.0)*self.rho[0,:]*self.u[0,0,:] -
                    (1.0/2.0)*self.rho[0,:]*self.u[1,0,:] )

    ### ************************************************
    ### Zou-He right wall velocity b.c.
    def zou_he_right_wall_velocity(self):

        lx = self.lx
        ly = self.ly

        self.u[0,lx,:] = self.u_right[0,:]
        self.u[1,lx,:] = self.u_right[1,:]

        self.rho[lx,:] = (self.g[0,lx,:] +
                          self.g[3,lx,:] +
                          self.g[4,lx,:] +
                      2.0*self.g[1,lx,:] +
                      2.0*self.g[5,lx,:] +
                      2.0*self.g[8,lx,:])/(1.0 + self.u[0,lx,:])

        self.g[2,lx,:] = (self.g[1,lx,:] -
                          (2.0/3.0)*self.rho[lx,:]*self.u[0,lx,:])

        self.g[6,lx,:] = (     self.g[5,lx,:]  +
                          0.5*(self.g[3,lx,:]  -
                               self.g[4,lx,:]) -
                          (1.0/6.0)*self.rho[lx,:]*self.u[0,lx,:] -
                          (1.0/2.0)*self.rho[lx,:]*self.u[1,lx,:] )

        self.g[7,lx,:] = (     self.g[8,lx,:]  -
                          0.5*(self.g[3,lx,:]  -
                               self.g[4,lx,:]) -
                          (1.0/6.0)*self.rho[lx,:]*self.u[0,lx,:] +
                          (1.0/2.0)*self.rho[lx,:]*self.u[1,lx,:] )

    ### ************************************************
    ### Zou-He right wall pressure b.c.
    def zou_he_right_wall_pressure(self):

        lx = self.lx
        ly = self.ly

        self.rho[lx,:] = self.rho_right[:]
        self.u[1,lx,:] = self.u_right[1,:]

        self.u[0,lx,:] =    (self.g[0,lx,:] +
                             self.g[3,lx,:] +
                             self.g[4,lx,:] +
                         2.0*self.g[1,lx,:] +
                         2.0*self.g[5,lx,:] +
                         2.0*self.g[8,lx,:])/self.rho[lx,:] - 1.0

        self.g[2,lx,:] = (self.g[1,lx,:] -
                          (2.0/3.0)*self.rho[lx,:]*self.u[0,lx,:])

        self.g[6,lx,:] = (     self.g[5,lx,:]  +
                          0.5*(self.g[3,lx,:]  -
                               self.g[4,lx,:]) -
                          (1.0/6.0)*self.rho[lx,:]*self.u[0,lx,:] -
                          (1.0/2.0)*self.rho[lx,:]*self.u[1,lx,:] )

        self.g[7,lx,:] = (     self.g[8,lx,:]  -
                          0.5*(self.g[3,lx,:]  -
                               self.g[4,lx,:]) -
                          (1.0/6.0)*self.rho[lx,:]*self.u[0,lx,:] +
                          (1.0/2.0)*self.rho[lx,:]*self.u[1,lx,:] )

    ### ************************************************
    ### Zou-He no-slip top wall velocity b.c.
    def zou_he_top_wall_velocity(self):

        lx = self.lx
        ly = self.ly

        self.u[0,:,ly] = self.u_top[0,:]
        self.u[1,:,ly] = self.u_top[1,:]

        self.rho[:,0] = (self.g[0,:,0] +
                         self.g[1,:,0] +
                         self.g[2,:,0] +
                     2.0*self.g[3,:,0] +
                     2.0*self.g[5,:,0] +
                     2.0*self.g[7,:,0] )/(1.0 + self.u[1,:,ly])

        self.g[4,:,ly] = (self.g[3,:,ly] -
                          (2.0/3.0)*self.rho[:,ly]*self.u[1,:,ly])

        self.g[8,:,ly] = (     self.g[7,:,ly]  -
                          0.5*(self.g[1,:,ly]  -
                               self.g[2,:,ly]) +
                          (1.0/2.0)*self.rho[:,ly]*self.u[0,:,ly] -
                          (1.0/6.0)*self.rho[:,ly]*self.u[1,:,ly] )

        self.g[6,:,ly] = (     self.g[5,:,ly]  +
                          0.5*(self.g[1,:,ly]  -
                               self.g[2,:,ly]) -
                          (1.0/2.0)*self.rho[:,ly]*self.u[0,:,ly] -
                          (1.0/6.0)*self.rho[:,ly]*self.u[1,:,ly] )

    ### ************************************************
    ### Zou-He no-slip bottom wall velocity b.c.
    def zou_he_bottom_wall_velocity(self):

        lx = self.lx
        ly = self.ly

        self.u[0,:,0] = self.u_bot[0,:]
        self.u[1,:,0] = self.u_bot[1,:]

        self.rho[:,0] = (self.g[0,:,0] +
                         self.g[1,:,0] +
                         self.g[2,:,0] +
                     2.0*self.g[4,:,0] +
                     2.0*self.g[6,:,0] +
                     2.0*self.g[8,:,0] )/(1.0 - self.u[1,:,0])

        self.g[3,:,0] = (self.g[4,:,0] +
                         (2.0/3.0)*self.rho[:,0]*self.u[1,:,0])

        self.g[5,:,0] = (     self.g[6,:,0]  -
                         0.5*(self.g[1,:,0]  -
                              self.g[2,:,0]) +
                         (1.0/2.0)*self.rho[:,0]*self.u[0,:,0] +
                         (1.0/6.0)*self.rho[:,0]*self.u[1,:,0] )

        self.g[7,:,0] = (     self.g[8,:,0]  +
                         0.5*(self.g[1,:,0]  -
                              self.g[2,:,0]) -
                         (1.0/2.0)*self.rho[:,0]*self.u[0,:,0] +
                         (1.0/6.0)*self.rho[:,0]*self.u[1,:,0] )

    ### ************************************************
    ### Zou-He bottom left corner
    def zou_he_bottom_left_corner(self):

        lx = self.lx
        ly = self.ly

        self.u[0,0,0] = self.u[0,1,0]
        self.u[1,0,0] = self.u[1,1,0]
        self.rho[0,0] = self.rho[1,0]

        self.g[1,0,0] = (self.g[2,0,0] +
                         (2.0/3.0)*self.rho[0,0]*self.u[0,0,0])

        self.g[3,0,0] = (self.g[4,0,0] +
                         (2.0/3.0)*self.rho[0,0]*self.u[1,0,0])

        self.g[5,0,0] = (self.g[6,0,0] +
                         (1.0/6.0)*self.rho[0,0]*self.u[0,0,0] +
                         (1.0/6.0)*self.rho[0,0]*self.u[1,0,0] )

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

        lx = self.lx
        ly = self.ly

        self.u[0,0,ly] = self.u[0,1,ly]
        self.u[1,0,ly] = self.u[1,1,ly]
        self.rho[0,ly] = self.rho[1,ly]

        self.g[1,0,ly] = (self.g[2,0,ly] +
                         (2.0/3.0)*self.rho[0,ly]*self.u[0,0,ly])

        self.g[4,0,ly] = (self.g[3,0,ly] -
                         (2.0/3.0)*self.rho[0,ly]*self.u[1,0,ly])

        self.g[8,0,ly] = (self.g[7,0,ly] +
                         (1.0/6.0)*self.rho[0,ly]*self.u[0,0,ly] -
                         (1.0/6.0)*self.rho[0,ly]*self.u[1,0,ly])


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

        self.u[0,lx,ly] = self.u[0,lx-1,ly]
        self.u[1,lx,ly] = self.u[1,lx-1,ly]
        self.rho[lx,ly] = self.rho[lx-1,ly]

        self.g[2,lx,ly] = (self.g[1,lx,ly] -
                         (2.0/3.0)*self.rho[lx,ly]*self.u[0,lx,ly])

        self.g[4,lx,ly] = (self.g[3,lx,ly] -
                         (2.0/3.0)*self.rho[lx,ly]*self.u[1,lx,ly])

        self.g[6,lx,ly] = (self.g[5,lx,ly] -
                         (1.0/6.0)*self.rho[lx,ly]*self.u[0,lx,ly] -
                         (1.0/6.0)*self.rho[lx,ly]*self.u[1,lx,ly])

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
        ly = self.ly

        self.u[0,lx,0] = self.u[0,lx-1,0]
        self.u[1,lx,0] = self.u[1,lx-1,0]
        self.rho[lx,0] = self.rho[lx-1,0]

        self.g[2,lx,0] = (self.g[1,lx,0] -
                         (2.0/3.0)*self.rho[lx,0]*self.u[0,lx,0])

        self.g[3,lx,0] = (self.g[4,lx,0] +
                         (2.0/3.0)*self.rho[lx,0]*self.u[1,lx,0])

        self.g[7,lx,0] = (self.g[8,lx,0] -
                         (1.0/6.0)*self.rho[lx,0]*self.u[0,lx,0] +
                         (1.0/6.0)*self.rho[lx,0]*self.u[1,lx,0])

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
    ### Output 2D flow amplitude
    def output_fields(self, it, freq, *args, **kwargs):

        # Handle inputs
        u_norm   = kwargs.get('u_norm',   True)
        u_stream = kwargs.get('u_stream', True)

        # Exit if no plotting
        if (it%freq != 0): return

        # Compute norm
        v = self.u.copy()
        v = np.sqrt(v[0,:,:]**2+v[1,:,:]**2)

        # Plot u norm
        if (u_norm):
            plt.clf()
            plt.imshow(np.rot90(v),
                       cmap = 'RdBu',
                       vmin = 0,
                       vmax = self.u_lbm,
                       interpolation = 'nearest')

            filename = self.png_dir+'u_norm_'+str(self.output_it)+'.png'
            plt.axis('off')
            plt.savefig(filename,
                        dpi=self.dpi)
            self.trim_white(filename)

        # Plot u streamlines
        if (u_stream):
            plt.clf()
            x = np.linspace(0, 1, self.nx)
            y = np.linspace(0, 1, self.ny)
            u = np.linspace(0, 1, 50)
            g = np.meshgrid(u,u)
            str_pts = list(zip(*(x.flat for x in g)))
            plt.streamplot(y, x,
                           self.u[1],
                           self.u[0],
                           linewidth = 0.2,
                           color='k',
                           arrowstyle = '-',
                           start_points=str_pts,
                           density = 50)

            filename = self.output_dir+'u_stream'+'.png'
            plt.axis('off')
            plt.gca().set_aspect('equal')
            plt.savefig(filename,
                        dpi=self.dpi)
            self.trim_white(filename)

        # Update counter
        self.output_it += 1

    ### ************************************************
    ### Add obstacle
    def add_obstacle(self, polygon):

        # Compute polygon bnds
        poly_bnds    = np.zeros((4))
        poly_bnds[0] = np.amin(polygon[:,0])
        poly_bnds[1] = np.amax(polygon[:,0])
        poly_bnds[2] = np.amin(polygon[:,1])
        poly_bnds[3] = np.amax(polygon[:,1])

        # Declare lattice arrays
        obstacle = np.empty((0,2), dtype=int)
        boundary = np.empty((0,3), dtype=int)
        ibb      = np.empty((1),   dtype=float)

        # Fill lattice
        print('Generating lattice...')
        for i in range(self.nx):
            for j in range(self.ny):
                pt = self.lattice_coords(i, j)

                # Check if pt is inside polygon bbox
                if ((pt[0] > poly_bnds[0]) and
                    (pt[0] < poly_bnds[1]) and
                    (pt[1] > poly_bnds[2]) and
                    (pt[1] < poly_bnds[3])):

                    if (self.is_inside(polygon, pt)):
                        self.lattice[i,j] = True
                        obstacle = np.append(obstacle,
                                             np.array([[i,j]]), axis=0)

        print('Found '+str(obstacle.shape[0])+' locations in obstacle')

        # Check area of obstacle
        area = 0.0

        for i in range(self.nx):
            for j in range(self.ny):
                if (self.lattice[i,j]): area += self.dx**2

        print('Obstacle area: '+str(area))

        # Build boundary of obstacle, i.e. 1st layer of fluid
        for k in range(len(obstacle)):
            i = obstacle[k,0]
            j = obstacle[k,1]

            for q in range(1,9):
                qb  = self.ns[q]
                cx  = self.c[q,0]
                cy  = self.c[q,1]
                ii  = i + cx
                jj  = j + cy

                if (not self.lattice[ii,jj]):
                    boundary = np.append(boundary,
                                         np.array([[ii,jj,qb]]), axis=0)

        # Some cells were counted multiple times, unique-sort them
        boundary = np.unique(boundary, axis=0)

        # Compute lattice-boundary distances if IBB is True
        if (self.IBB):
            for k in range(len(boundary)):
                i    = boundary[k,0]
                j    = boundary[k,1]
                q    = boundary[k,2]
                pt   = self.lattice_coords(i, j)
                x    = polygon[:,0] - pt[0]
                y    = polygon[:,1] - pt[1]
                dist = np.sqrt(np.square(x) + np.square(y))
                mpt  = np.argmin(dist)
                mdst = dist[mpt]/(self.dx*np.linalg.norm(self.c[q]))
                ibb  = np.append(ibb, mdst)

        # Add obstacle
        obs = Obstacle(polygon, area, boundary, ibb)
        self.obstacles.append(obs)

        # Printings
        print('Found '+str(boundary.shape[0])+' on obstacle boundary')

    ### ************************************************
    ### Get lattice coordinates from integers
    def lattice_coords(self, i, j):

        # Compute and return the coordinates of the lattice node (i,j)
        dx = (self.x_max - self.x_min)/(self.nx - 1)
        dy = (self.y_max - self.y_min)/(self.ny - 1)
        x  = self.x_min + i*dx
        y  = self.y_min + j*dy

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

        for obs in range(len(self.obstacles)):
            for k in range(len(self.obstacles[obs].boundary)):
                i = self.obstacles[obs].boundary[k,0]
                j = self.obstacles[obs].boundary[k,1]
                lat[i,j] += 2.0

        # Plot and save image of lattice
        filename = self.output_dir+self.name+'.png'

        plt.imsave(filename,
                   np.rot90(lat),
                   vmin=-1.0,
                   vmax= 2.0)
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
    ### Set inlet poiseuille fields
    def set_inlet_poiseuille(self, u_lbm, rho_lbm, it, sigma):

        lx              = self.lx
        ly              = self.ly

        self.u_left     = np.zeros((2, self.ny))
        self.u_right    = np.zeros((2, self.ny))
        self.u_top      = np.zeros((2, self.nx))
        self.u_bot      = np.zeros((2, self.nx))
        self.rho_right  = np.ones(self.ny)
        self.rho_right *= rho_lbm

        for j in range(self.ny):
            pt               = self.lattice_coords(0, j)
            self.u_left[:,j] = u_lbm*self.poiseuille(pt, it, sigma)

    ### ************************************************
    ### Set full poiseuille fields
    def set_full_poiseuille(self, u_lbm, rho_lbm):

        lx              = self.lx
        ly              = self.ly

        self.u_left     = np.zeros((2, self.ny))
        self.u_right    = np.zeros((2, self.ny))
        self.u_top      = np.zeros((2, self.nx))
        self.u_bot      = np.zeros((2, self.nx))
        self.rho_right  = np.ones(self.ny)
        self.rho_right *= rho_lbm

        for j in range(self.ny):
            for i in range(self.nx):
                pt               = self.lattice_coords(i, j)
                u                = u_lbm*self.poiseuille(pt, 1, 1.0e-10)
                self.u_left[:,j] = u
                self.u[:,i,j]    = u

    ### ************************************************
    ### Set driven cavity fields
    def set_cavity(self, u_lbm):

        lx              = self.lx
        ly              = self.ly

        self.u_left     = np.zeros((2, self.ny))
        self.u_right    = np.zeros((2, self.ny))
        self.u_top      = np.zeros((2, self.nx))
        self.u_bot      = np.zeros((2, self.nx))

        self.u_top[0,:] = u_lbm
        self.u[0,:,ly]  = self.u_top[0,:]
        self.u[1,:,ly]  = self.u_top[1,:]

    ### ************************************************
    ### Poiseuille flow
    def poiseuille(self, pt, it, sigma):

        x    = pt[0]
        y    = pt[1]
        H    = self.y_max - self.y_min
        u    = np.zeros(2)
        u[0] = 4.0*(self.y_max-y)*(y-self.y_min)/H**2

        val  = it
        ret  = (1.0 - math.exp(-val**2/(2.0*sigma**2)))
        u   *= ret

        return u

    ### ************************************************
    ### Poiseuille error in the middle of the domain
    def poiseuille_error(self, u_lbm):

        u_error = np.zeros((2,self.ny))
        nx      = math.floor(self.nx/2)

        for j in range(self.ny):
            pt   = self.lattice_coords(nx,j)
            u_ex = self.poiseuille(pt, 1.0e10, 1)
            u    = self.u[:,nx,j]

            u_error[0,j] = u[0]/u_lbm
            u_error[1,j] = u_ex[0]

        # Write to file
        filename = self.output_dir+'poiseuille'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {} {}\n'.format(j*self.dx,
                                            u_error[0,j],
                                            u_error[1,j]))

    ### ************************************************
    ### Cavity error in the middle of the domain
    def cavity_error(self, u_lbm):

        ux_error = np.zeros((self.nx))
        uy_error = np.zeros((self.ny))
        nx       = math.floor(self.nx/2)
        ny       = math.floor(self.ny/2)

        for i in range(self.nx):
            uy_error[i] = self.u[1,i,ny]/u_lbm

        for j in range(self.ny):
            ux_error[j] = self.u[0,nx,j]/u_lbm

        # Write to files
        filename = self.output_dir+'cavity_uy'
        with open(filename, 'w') as f:
            for i in range(self.nx):
                f.write('{} {}\n'.format(i*self.dx, uy_error[i]))
        filename = self.output_dir+'cavity_ux'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {}\n'.format(j*self.dx, ux_error[j]))

    ### ************************************************
    ### Check stopping criterion
    def check_stop(self):

        if (self.stop == 'it'):
            if (self.it > self.it_max):
                self.compute = False
                print('\n')

        if (self.stop == 'obs'):
            if (self.drag_buff.obs_cv and self.lift_buff.obs_cv):
                self.compute = False
                print('\n')

        self.it += 1

    ### ************************************************
    ### Iteration printings
    def it_printings(self):

        if (self.stop == 'it'):
            print('# it = '+str(self.it)+' / '+str(self.it_max),
                  end='\r')
        if (self.stop == 'obs'):
            str_d = "{:10.6f}".format(self.drag_buff.obs)
            str_l = "{:10.6f}".format(self.lift_buff.obs)

            print('it = '+str(self.it)+
                  ', avg drag = '+str_d+', avg lift = '+str_l, end='\r')

### ************************************************
### Compute equilibrium state
@jit(nopython=True,parallel=True,cache=True)
def nb_equilibrium(u, c, w, rho, g_eq):

    # Compute velocity term
    v = (3.0/2.0)*(u[0,:,:]**2 + u[1,:,:]**2)

    # Compute equilibrium
    for q in range(9):
        t                 = 3.0*(u[0,:,:]*c[q,0] +
                                 u[1,:,:]*c[q,1])
        g_eq[q,:,:]  = (1.0 + t + 0.5*t**2 - v)
        g_eq[q,:,:] *= rho[:,:]*w[q]

### ************************************************
### Collision and streaming
@jit(nopython=True,parallel=True,cache=True)
def nb_col_str(g, g_eq, g_up, om_p, om_m, c, ns, nx, ny, lx, ly):

    # Take care of q=0 first
    g_up[0,:,:] = g[0,:,:] - om_p*(g[0,:,:] - g_eq[0,:,:])
    g   [0,:,:] = g_up[0,:,:]

    # Collide other indices
    for q in range(1,9):
        qb = ns[q]

        g_up[q,:,:] = (g   [q,:,:]
                       - om_p*0.5*(g   [q,:,:]
                                   + g   [qb,:,:]
                                   - g_eq[q,:,:]
                                   - g_eq[qb,:,:])
                       - om_m*0.5*(g   [q,:,:]
                                   - g   [qb,:,:]
                                   - g_eq[q,:,:]
                                   + g_eq[qb,:,:]))

    # Stream
    g[1,1:nx, :  ] = g_up[1,0:lx, :  ]
    g[2,0:lx, :  ] = g_up[2,1:nx, :  ]
    g[3, :,  1:ny] = g_up[3, :,  0:ly]
    g[4, :,  0:ly] = g_up[4, :,  1:ny]
    g[5,1:nx,1:ny] = g_up[5,0:lx,0:ly]
    g[6,0:lx,0:ly] = g_up[6,1:nx,1:ny]
    g[7,0:lx,1:ny] = g_up[7,1:nx,0:ly]
    g[8,1:nx,0:ly] = g_up[8,0:lx,1:ny]
