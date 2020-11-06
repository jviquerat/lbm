# Generic imports
import math
import time

# Custom imports
from shapes  import *
from lattice import *

###############################################
### LBM poiseuille with array of obstacles
###############################################

# Obstacles
obs_name     = 'main'
obs_npts     = 4
obs_nspts    = 100
obs_type     = 'square'
obs_size     = 0.1
radius       = 0.5

n_obs        = 8
obs_pos      = []
for i in range(n_obs):
    pos = [radius*math.cos(2.0*math.pi*float(i)/n_obs),
           radius*math.sin(2.0*math.pi*float(i)/n_obs)]
    obs_pos.append(pos)

# Domain size
x_min       =-1.0
x_max       = 8.0
y_min       =-1.0
y_max       = 1.0

# Free parameters
# L_lbm correponds to y length
# u_lbm corresponds to max velocity
Re_lbm      = 100.0
u_lbm       = 0.03
L_lbm       = 400

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
ny          = L_lbm
u_avg       = 2.0*u_lbm/3.0
D_lbm       = math.floor(ny*obs_size/(y_max-y_min))
nu_lbm      = u_avg*D_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
dx          = (y_max-y_min)/ny
dy          = dx
nx          = math.floor(ny*(x_max-x_min)/(y_max-y_min))

# Other parameters
output_freq = 1000
dpi         = 300
IBB         = True
stop        = 'obs'
obs_cv_ct   = 1.0e-2
obs_cv_nb   = 1000
t_max       = 0.02
it_max      = math.floor(t_max/dt)
sigma       = math.floor(10*nx)

# Initialize lattice
lat = Lattice(nx        = nx,
              ny        = ny,
              dx        = dx,
              dt        = dt,
              tau_lbm   = tau_lbm,
              Re_lbm    = Re_lbm,
              u_lbm     = u_lbm,
              L_lbm     = D_lbm,
              nu_lbm    = nu_lbm,
              rho_lbm   = rho_lbm,
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

# Generate obstacles
obs = [None]*n_obs
for i in range(n_obs):
    obs[i] = generate_shape(obs_npts,
                            obs_pos[i],
                            obs_type,
                            obs_size,
                            obs_name,
                            obs_nspts,
                            lat.output_dir)
    lat.add_obstacle(obs[i].curve_pts, i+1)
lat.generate_image()

# Initialize fields
lat.set_inlet_poiseuille(u_lbm, rho_lbm, 0, sigma)
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

    # Progressively impose Poiseuille
    lat.set_inlet_poiseuille(u_lbm, rho_lbm, lat.it, sigma)

    # Compute macroscopic fields
    lat.macro()

    # Output field
    lat.output_fields(lat.it,
                      output_freq,
                      u_norm   = True,
                      u_stream = False)

    # Compute equilibrium state
    lat.equilibrium()

    # Collisions
    lat.collision_stream()

    # Boundary conditions
    for i in range(n_obs):
        lat.bounce_back_obstacle(i)
    lat.zou_he_bottom_wall_velocity()
    lat.zou_he_left_wall_velocity()
    lat.zou_he_top_wall_velocity()
    lat.zou_he_right_wall_pressure()
    lat.zou_he_bottom_left_corner()
    lat.zou_he_top_left_corner()
    lat.zou_he_top_right_corner()
    lat.zou_he_bottom_right_corner()

    # Compute drag/lift
    drag = 0.0
    lift = 0.0
    for i in range(n_obs):
        drag_obs, lift_obs = lat.drag_lift(i, rho_lbm, u_avg, D_lbm)
        drag              += drag_obs
        lift              += lift_obs
    lat.add_buff(drag, lift, lat.it)

    # Check stopping criterion
    lat.check_stop()

# Count time
end_time = time.time()
print("# Loop time = {:f}".format(end_time - start_time))
