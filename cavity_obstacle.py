# Generic imports
import math
import time

# Custom imports
from shapes  import *
from lattice import *

###############################################
### LBM lid-driven cavity
###############################################

### Free parameters
Re_lbm      = 1000.0
u_lbm       = 0.03
L_lbm       = 100

# Obstacle parameters
pos  = [0.7, 0.7]
size = 0.1

# Domain size
x_min       =-1.0
x_max       = 1.0
y_min       =-1.0
y_max       = 1.0

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
ny          = L_lbm
D_lbm       = math.floor(ny*size/(y_max-y_min))
nu_lbm      = u_lbm*L_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
dx          = (y_max-y_min)/ny
dy          = dx
nx          = math.floor(ny*(x_max-x_min)/(y_max-y_min))

# Other parameters
output_freq = 2000
dpi         = 200
IBB         = True
stop        = 'obs'
obs_cv_ct   = 1.0e-3
obs_cv_nb   = 5000
t_max       = 0.02
it_max      = math.floor(t_max/dt)
sigma       = math.floor(5*nx)

# Initialize lattice
lat = Lattice(nx      = nx,
              ny      = nx,
              dx      = 1.0/nx,
              dt      = dt,
              tau_lbm = tau_lbm,
              Re_lbm  = Re_lbm,
              u_lbm   = u_lbm,
              L_lbm   = L_lbm,
              nu_lbm  = nu_lbm,
              rho_lbm = rho_lbm,
              x_min     = x_min,
              x_max     = x_max,
              y_min     = y_min,
              y_max     = y_max,
              dpi       = dpi,
              IBB       = IBB,
              stop      = stop,
              t_max     = t_max,
              it_max    = it_max,
              obs_cv_ct = obs_cv_ct,
              obs_cv_nb = obs_cv_nb)

# Generate obstacle
obs = generate_shape(4, pos, 'square', size, 'main', 100, lat.output_dir)
lat.add_obstacle(obs.curve_pts, 1)
lat.generate_image()

# Initialize fields
lat.rho *= rho_lbm
lat.u[:,np.where(lat.lattice > 0.0)] = 0.0

# Set initial distributions
lat.equilibrium()
lat.g = lat.g_eq.copy()

# Count time
start_time = time.time()

# Solve
print('### Solving')
while (lat.compute):

    # Printings
    lat.it_printings()

    # Impose velocity
    lat.set_cavity(u_lbm, -u_lbm, 0.0, 0.0, lat.it, sigma)

    # Compute macroscopic fields
    lat.macro()

    # Output field
    lat.output_fields(lat.it,
                          output_freq,
                          u_norm   = True,
                          u_stream = False)

    # Compute equilibrium state
    lat.equilibrium()

    # Streaming
    lat.collision_stream()

    # Boundary conditions
    lat.bounce_back_obstacle(0)
    lat.zou_he_bottom_wall_velocity()
    lat.zou_he_left_wall_velocity()
    lat.zou_he_right_wall_velocity()
    lat.zou_he_top_wall_velocity()
    lat.zou_he_bottom_left_corner()
    lat.zou_he_top_left_corner()
    lat.zou_he_top_right_corner()
    lat.zou_he_bottom_right_corner()

    # Compute drag/lift
    drag, lift = lat.drag_lift(0, rho_lbm, u_lbm, D_lbm)
    lat.add_buff(drag, lift, lat.it)

    # Check stopping criterion
    lat.check_stop()

# Count time
end_time = time.time()
print("# Loop time = {:f}".format(end_time - start_time))

# Output streamlines
lat.output_fields(1, 1,
                  u_norm   = False,
                  u_stream = True)
