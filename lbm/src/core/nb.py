# Generic imports
import numba as     nb
from   numba import jit, njit

### ************************************************
### Compute equilibrium state
@njit(parallel=True,cache=True,fastmath=True,nogil=True)
def nb_equilibrium(u, c, w, rho, g_eq):

    # Compute velocity term
    v = 1.5*(u[0,:,:]**2 + u[1,:,:]**2)

    # Compute equilibrium
    for q in nb.prange(9):
        t            = 3.0*(u[0,:,:]*c[q,0] + u[1,:,:]*c[q,1])
        g_eq[q,:,:]  = (1.0 + t + 0.5*t**2 - v)
        g_eq[q,:,:] *= rho[:,:]*w[q]

### ************************************************
### Collision and streaming
@njit(parallel=True,cache=True,fastmath=True,nogil=True)
def nb_col_str(g, g_eq, g_up, om_p, om_m, c, ns, nx, ny, lx, ly):

    # Take care of q=0 first
    g_up[0,:,:] = (1.0-om_p)*g[0,:,:] + om_p*g_eq[0,:,:]
    g   [0,:,:] = g_up[0,:,:]

    # Collide other indices
    for q in nb.prange(1,9):
        qb = ns[q]

        g_up[q,:,:] = ((1.0-0.5*(om_p+om_m))*g[q,:,:]   -
                            0.5*(om_p-om_m)*g[qb,:,:]   +
                            0.5*(om_p+om_m)*g_eq[q,:,:] +
                            0.5*(om_p-om_m)*g_eq[qb,:,:])

    # Stream
    g[1,1:nx, :  ] = g_up[1,0:lx, :  ]
    g[2,0:lx, :  ] = g_up[2,1:nx, :  ]
    g[3, :,  1:ny] = g_up[3, :,  0:ly]
    g[4, :,  0:ly] = g_up[4, :,  1:ny]
    g[5,1:nx,1:ny] = g_up[5,0:lx,0:ly]
    g[6,0:lx,0:ly] = g_up[6,1:nx,1:ny]
    g[7,0:lx,1:ny] = g_up[7,1:nx,0:ly]
    g[8,1:nx,0:ly] = g_up[8,0:lx,1:ny]

### ************************************************
### Compute drag and lift
@njit(parallel=True,cache=True,fastmath=True,nogil=True)
def nb_drag_lift(boundary, ns, c, g_up, g, R_ref, U_ref, L_ref):

    # Initialize
    fx     = 0.0
    fy     = 0.0

    # Loop over obstacle array
    for k in nb.prange(len(boundary)):
        i   = boundary[k,0]
        j   = boundary[k,1]
        q   = boundary[k,2]
        qb  = ns[q]
        cx  = c[q,0]
        cy  = c[q,1]
        g0  = g_up[q,i,j] + g[qb,i,j]

        fx += g0*cx
        fy += g0*cy

    # Normalize coefficient
    Cx =-2.0*fx/(R_ref*L_ref*U_ref**2)
    Cy =-2.0*fy/(R_ref*L_ref*U_ref**2)

    return Cx, Cy

### ************************************************
### Obstacle halfway bounce-back no-slip b.c.
@njit(parallel=True,cache=True,nogil=True)
def nb_bounce_back_obstacle(IBB, boundary, ns, sc,
                            obs_ibb, g_up, g, u, lattice):

    # Interpolated BB
    if (IBB):
        for k in nb.prange(len(boundary)):
            i  = boundary[k,0]
            j  = boundary[k,1]
            q  = boundary[k,2]
            qb = ns[q]
            c  = sc[q,:]
            cb = sc[qb,:]
            im = i + cb[0]
            jm = j + cb[1]
            imm = i + 2*cb[0]
            jmm = j + 2*cb[1]

            p  = obs_ibb[k]
            pp = 2.0*p
            if (p < 0.5):
                g[qb,i,j] = (p*(pp+1.0)*g_up[q,i,j]
                             + (1.0+pp)*(1.0-pp)*g_up[q,im,jm]
                             - p*(1.0-pp)*g_up[q,imm,jmm])
            else:
                g[qb,i,j] = ((1.0/(p*(pp+1.0)))*g_up[q,i,j] +
                             ((pp-1.0)/p)*g_up[qb,i,j] +
                             ((1.0-pp)/(1.0+pp))*g_up[qb,im,jm])

    # Regular BB
    if (not IBB):
        for k in nb.prange(len(boundary)):
            i  = boundary[k,0]
            j  = boundary[k,1]
            q  = boundary[k,2]
            qb = ns[q]
            c  = sc[q,:]
            ii = i + c[0]
            jj = j + c[1]

            g[qb,i,j] = g_up[q,i,j]

### ************************************************
### Zou-He left wall velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_left_wall_velocity(lx, ly, u, u_left, rho, g):

    cst1 = 2.0/3.0
    cst2 = 1.0/6.0
    cst3 = 1.0/2.0

    u[0,0,:] = u_left[0,:]
    u[1,0,:] = u_left[1,:]

    rho[0,:] = (g[0,0,:] + g[3,0,:] + g[4,0,:] +
                2.0*g[2,0,:] + 2.0*g[6,0,:] +
                2.0*g[7,0,:] )/(1.0 - u[0,0,:])

    g[1,0,:] = (g[2,0,:] + cst1*rho[0,:]*u[0,0,:])

    g[5,0,:] = (g[6,0,:] - cst3*(g[3,0,:] - g[4,0,:]) +
                cst2*rho[0,:]*u[0,0,:] +
                cst3*rho[0,:]*u[1,0,:] )

    g[8,0,:] = (g[7,0,:] + cst3*(g[3,0,:] - g[4,0,:]) +
                cst2*rho[0,:]*u[0,0,:] -
                cst3*rho[0,:]*u[1,0,:] )

### ************************************************
### Zou-He right wall velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_right_wall_velocity(lx, ly, u, u_right, rho, g):

    cst1 = 2.0/3.0
    cst2 = 1.0/6.0
    cst3 = 1.0/2.0

    u[0,lx,:] = u_right[0,:]
    u[1,lx,:] = u_right[1,:]

    rho[lx,:] = (g[0,lx,:] + g[3,lx,:] + g[4,lx,:] +
                 2.0*g[1,lx,:] + 2.0*g[5,lx,:] +
                 2.0*g[8,lx,:])/(1.0 + u[0,lx,:])

    g[2,lx,:] = (g[1,lx,:] - cst1*rho[lx,:]*u[0,lx,:])

    g[6,lx,:] = (g[5,lx,:] + cst3*(g[3,lx,:] - g[4,lx,:]) -
                 cst2*rho[lx,:]*u[0,lx,:] -
                 cst3*rho[lx,:]*u[1,lx,:] )

    g[7,lx,:] = (g[8,lx,:] - cst3*(g[3,lx,:] - g[4,lx,:]) -
                 cst2*rho[lx,:]*u[0,lx,:] +
                 cst3*rho[lx,:]*u[1,lx,:] )

### ************************************************
### Zou-He right wall pressure b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_right_wall_pressure(lx, ly, u, rho_right, u_right, rho, g):

    cst1 = 2.0/3.0
    cst2 = 1.0/6.0
    cst3 = 1.0/2.0

    rho[lx,:] = rho_right[:]
    u[1,lx,:] = u_right[1,:]

    u[0,lx,:] = (g[0,lx,:] + g[3,lx,:] + g[4,lx,:] +
                 2.0*g[1,lx,:] + 2.0*g[5,lx,:] +
                 2.0*g[8,lx,:])/rho[lx,:] - 1.0

    g[2,lx,:] = (g[1,lx,:] - cst1*rho[lx,:]*u[0,lx,:])

    g[6,lx,:] = (g[5,lx,:] + cst3*(g[3,lx,:] - g[4,lx,:]) -
                 cst2*rho[lx,:]*u[0,lx,:] -
                 cst3*rho[lx,:]*u[1,lx,:] )

    g[7,lx,:] = (g[8,lx,:] - cst3*(g[3,lx,:] - g[4,lx,:]) -
                 cst2*rho[lx,:]*u[0,lx,:] +
                 cst3*rho[lx,:]*u[1,lx,:] )

### ************************************************
### Zou-He no-slip top wall velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_top_wall_velocity(lx, ly, u, u_top, rho, g):

    cst1 = 2.0/3.0
    cst2 = 1.0/6.0
    cst3 = 1.0/2.0

    u[0,:,ly] = u_top[0,:]
    u[1,:,ly] = u_top[1,:]

    rho[:,ly] = (g[0,:,ly] + g[1,:,ly] + g[2,:,ly] +
                2.0*g[3,:,ly] + 2.0*g[5,:,ly] +
                2.0*g[7,:,ly])/(1.0 + u[1,:,ly])

    g[4,:,ly] = (g[3,:,ly] - cst1*rho[:,ly]*u[1,:,ly])

    g[8,:,ly] = (g[7,:,ly] - cst3*(g[1,:,ly] - g[2,:,ly]) +
                 cst3*rho[:,ly]*u[0,:,ly] -
                 cst2*rho[:,ly]*u[1,:,ly] )

    g[6,:,ly] = (g[5,:,ly] + cst3*(g[1,:,ly] - g[2,:,ly]) -
                 cst3*rho[:,ly]*u[0,:,ly] -
                 cst2*rho[:,ly]*u[1,:,ly] )

### ************************************************
### Zou-He no-slip bottom wall velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_bottom_wall_velocity(lx, ly, u, u_bot, rho, g):

    cst1 = 2.0/3.0
    cst2 = 1.0/6.0
    cst3 = 1.0/2.0

    u[0,:,0] = u_bot[0,:]
    u[1,:,0] = u_bot[1,:]

    rho[:,0] = (g[0,:,0] + g[1,:,0] + g[2,:,0] +
                2.0*g[4,:,0] + 2.0*g[6,:,0] +
                2.0*g[8,:,0] )/(1.0 - u[1,:,0])

    g[3,:,0] = (g[4,:,0] + cst1*rho[:,0]*u[1,:,0])

    g[5,:,0] = (g[6,:,0] - cst3*(g[1,:,0] - g[2,:,0]) +
                cst3*rho[:,0]*u[0,:,0] +
                cst2*rho[:,0]*u[1,:,0] )

    g[7,:,0] = (g[8,:,0] + cst3*(g[1,:,0] - g[2,:,0]) -
                cst3*rho[:,0]*u[0,:,0] +
                cst2*rho[:,0]*u[1,:,0] )

### ************************************************
### Zou-He no-slip bottom left corner velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_bottom_left_corner_velocity(lx, ly, u, rho, g):

    u[0,0,0] = u[0,1,0]
    u[1,0,0] = u[1,1,0]

    rho[0,0] = rho[1,0]

    g[1,0,0] = (g[2,0,0] + (2.0/3.0)*rho[0,0]*u[0,0,0])

    g[3,0,0] = (g[4,0,0] + (2.0/3.0)*rho[0,0]*u[1,0,0])

    g[5,0,0] = (g[6,0,0] + (1.0/6.0)*rho[0,0]*u[0,0,0]
                         + (1.0/6.0)*rho[0,0]*u[1,0,0] )

    g[7,0,0] = 0.0
    g[8,0,0] = 0.0

    g[0,0,0] = (rho[0,0]
                - g[1,0,0] - g[2,0,0] - g[3,0,0] - g[4,0,0]
                - g[5,0,0] - g[6,0,0] - g[7,0,0] - g[8,0,0] )

### ************************************************
### Zou-He no-slip top left corner velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_top_left_corner_velocity(lx, ly, u, rho, g):

    u[0,0,ly] = u[0,1,ly]
    u[1,0,ly] = u[1,1,ly]

    rho[0,ly] = rho[1,ly]

    g[1,0,ly] = (g[2,0,ly] + (2.0/3.0)*rho[0,ly]*u[0,0,ly])

    g[4,0,ly] = (g[3,0,ly] - (2.0/3.0)*rho[0,ly]*u[1,0,ly])

    g[8,0,ly] = (g[7,0,ly] + (1.0/6.0)*rho[0,ly]*u[0,0,ly]
                           - (1.0/6.0)*rho[0,ly]*u[1,0,ly])


    g[5,0,ly] = 0.0
    g[6,0,ly] = 0.0

    g[0,0,ly] = (rho[0,ly]
                 - g[1,0,ly] - g[2,0,ly] - g[3,0,ly] - g[4,0,ly]
                 - g[5,0,ly] - g[6,0,ly] - g[7,0,ly] - g[8,0,ly] )

### ************************************************
### Zou-He no-slip top right corner velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_top_right_corner_velocity(lx, ly, u, rho, g):

    u[0,lx,ly] = u[0,lx-1,ly]
    u[1,lx,ly] = u[1,lx-1,ly]

    rho[lx,ly] = rho[lx-1,ly]

    g[2,lx,ly] = (g[1,lx,ly] - (2.0/3.0)*rho[lx,ly]*u[0,lx,ly])

    g[4,lx,ly] = (g[3,lx,ly] - (2.0/3.0)*rho[lx,ly]*u[1,lx,ly])

    g[6,lx,ly] = (g[5,lx,ly] - (1.0/6.0)*rho[lx,ly]*u[0,lx,ly]
                             - (1.0/6.0)*rho[lx,ly]*u[1,lx,ly])

    g[7,lx,ly] = 0.0
    g[8,lx,ly] = 0.0

    g[0,lx,ly] = (rho[lx,ly]
                  - g[1,lx,ly] - g[2,lx,ly] - g[3,lx,ly] - g[4,lx,ly]
                  - g[5,lx,ly] - g[6,lx,ly] - g[7,lx,ly] - g[8,lx,ly] )

### ************************************************
### Zou-He no-slip bottom right corner velocity b.c.
@njit(cache=True,fastmath=True,nogil=True)
def nb_zou_he_bottom_right_corner_velocity(lx, ly, u, rho, g):

    u[0,lx,0] = u[0,lx-1,0]
    u[1,lx,0] = u[1,lx-1,0]

    rho[lx,0] = rho[lx-1,0]

    g[2,lx,0] = (g[1,lx,0] - (2.0/3.0)*rho[lx,0]*u[0,lx,0])

    g[3,lx,0] = (g[4,lx,0] + (2.0/3.0)*rho[lx,0]*u[1,lx,0])

    g[7,lx,0] = (g[8,lx,0] - (1.0/6.0)*rho[lx,0]*u[0,lx,0]
                           + (1.0/6.0)*rho[lx,0]*u[1,lx,0])

    g[5,lx,0] = 0.0
    g[6,lx,0] = 0.0

    g[0,lx,0] = (rho[lx,0]
                 - g[1,lx,0] - g[2,lx,0] - g[3,lx,0] - g[4,lx,0]
                 - g[5,lx,0] - g[6,lx,0] - g[7,lx,0] - g[8,lx,0] )
