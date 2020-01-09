# Generic imports
import os
import datetime
from numpy        import *
from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib import cm

### General data
time        = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir = './results/'
output_dir  = results_dir+str(time)+'/'
png_dir     = output_dir+'./png/'

if (not os.path.exists(results_dir)): os.makedirs(results_dir)
if (not os.path.exists(output_dir)):  os.makedirs(output_dir)
if (not os.path.exists(png_dir)):     os.makedirs(png_dir)

### Flow parameters
maxIter = 20000     # total nb of iterations
Re      = 1000.0    # Re number
nx      = 520
ny      = 180
ly      = ny-1.0    # lattice dimension
q       = 9         # lattice population
cx      = nx/4      # cylinder coordinates
cy      = ny/2      # cylinder coordinates
r       = ny/9;     # cylinder radius
uLB     = 0.04      # velocity in lattice units
nulb    = uLB*r/Re
omega   = 1.0/(3.*nulb+0.5) # relaxation parameter.

### Lattice definition
# Set lattice velocities
c = array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # lattice velocities

# Set lattice weights
t                                      = ones(q);
t[asarray([norm(ci)<1.1 for ci in c])] = 1./9.  # cardinal       values
t[asarray([norm(ci)>1.1 for ci in c])] = 1./36. # extra-cardinal values
t[0]                                   = 4./9.  # central        value


noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)]
i1 = arange(q)[asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
i2 = arange(q)[asarray([ci[0]==0 for ci in c])] # Vertical middle.
i3 = arange(q)[asarray([ci[0]>0  for ci in c])] # Unknown on left wall.

print(noslip)

###### Function Definitions ####################################################
sumpop = lambda fin: sum(fin,axis=0) # Helper function for density computation.
def equilibrium(rho,u):              # Equilibrium distribution function.
    cu   = 3.0 * dot(c,u.transpose(1,0,2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = zeros((q,nx,ny))
    for i in range(q): feq[i,:,:] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq

###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
obstacle = fromfunction(lambda x,y: (x-cx)**2+(y-cy)**2<r**2, (nx,ny))
#vel      = fromfunction(lambda d,x,y: (1-d)*uLB*(1.0+1e-4*sin(y/ly*2*pi)),(2,nx,ny))
vel      = fromfunction(lambda d,x,y: (1-d)*uLB,(2,nx,ny))
feq      = equilibrium(1.0,vel)
fin      = feq.copy()

###### Main time loop ##########################################################
for time in range(maxIter):
    fin[i1,-1,:] = fin[i1,-2,:] # Right wall: outflow condition.
    rho          = sumpop(fin)  # Calculate macroscopic density and velocity.
    u            = dot(c.transpose(), fin.transpose((1,0,2)))/rho

    u[:,0,:]     = vel[:,0,:] # Left wall: compute density from known populations.
    rho[0,:]     = 1./(1.-u[0,0,:]) * (sumpop(fin[i2,0,:])+2.*sumpop(fin[i1,0,:]))

    feq = equilibrium(rho,u) # Left wall: Zou/He boundary condition.
    fin[i3,0,:] = fin[i1,0,:] + feq[i3,0,:] - fin[i1,0,:]
    fout = fin - omega * (fin - feq)  # Collision step.
    for i in range(q): fout[i,obstacle] = fin[noslip[i],obstacle]
    for i in range(q): # Streaming step.
        fin[i,:,:] = roll(roll(fout[i,:,:],c[i,0],axis=0),c[i,1],axis=1)

    if (time%100==0): # Visualization
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(),
                   cmap = cm.inferno,
                   vmin = 0.0,
                   vmax = 0.05)
        plt.savefig(png_dir+"vel."+str(time/100).zfill(4)+".png")
