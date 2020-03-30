# Generic imports
import math
import time

# Custom imports
from shapes_utils  import *
from lattice_utils import *

###############################################
### LBM poiseuille with obstacle
###############################################

# Count time
start_time = time.time()

# Shape1 parameters
shape1_name     = 'main'
shape1_npts     = 4
shape1_nspts    = 600
shape1_type     = 'square'
shape1_size     = 0.1
shape1_position = [0.0, 0.0]

# Domain size
x_min       =-0.2
x_max       = 2.0
y_min       =-0.2
y_max       = 0.22

# Free parameters
# L_lbm correponds to y length
# u_lbm corresponds to max velocity
Re_lbm      = 20.0
u_lbm       = 0.03
L_lbm       = 100

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
ny          = L_lbm
u_avg       = 2.0*u_lbm/3.0
D_lbm       = math.floor(ny*shape1_size/(y_max-y_min))
nu_lbm      = u_avg*D_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
dx          = (y_max-y_min)/ny
dy          = dx
nx          = math.floor(ny*(x_max-x_min)/(y_max-y_min))

# Other parameters
output_freq = 500
dpi         = 300
IBB         = True
stop        = 'it'
t_max       = 1.0
it_max      = math.floor(t_max/dt)
sigma       = math.floor(10*nx)

# Initialize lattice
lattice = Lattice(nx        = nx,
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
                  it_max    = it_max)

# Generate main shape and add to lattice
shape1 = generate_shape(shape1_npts,
                       shape1_position,
                       shape1_type,
                       shape1_size,
                       shape1_name,
                       shape1_nspts,
                       lattice.output_dir)
lattice.add_obstacle(shape1.curve_pts)
lattice.generate_image()

# Initialize fields
lattice.set_inlet_poiseuille(u_lbm, rho_lbm, 0, sigma)
lattice.u[:,lattice.lattice] = 0.0

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq.copy()

# Count time
end_time = time.time()
print("Pre-processing time = {}".format(end_time - start_time))
start_time = time.time()

# Solve
print('Solving...')
while (lattice.compute):

    lattice.it_printings()

    # Progressively impose Poiseuille
    lattice.set_inlet_poiseuille(u_lbm, rho_lbm, lattice.it, sigma)

    # Compute macroscopic fields
    lattice.macro()

    # Output field
    lattice.output_fields(lattice.it,
                          output_freq,
                          u_norm   = True,
                          u_stream = False)

    # Compute equilibrium state
    lattice.equilibrium()

    # Collisions
    lattice.collision_stream()

    # Boundary conditions
    lattice.bounce_back_obstacle(0)
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_top_wall_velocity()
    lattice.zou_he_right_wall_pressure()
    lattice.zou_he_bottom_left_corner()
    lattice.zou_he_top_left_corner()
    lattice.zou_he_top_right_corner()
    lattice.zou_he_bottom_right_corner()

    # Compute drag/lift
    lattice.drag_lift(0, lattice.it, rho_lbm, u_avg, D_lbm)

    # Check stopping criterion
    lattice.check_stop()

# Count time
end_time = time.time()
print("Loop time = {}".format(end_time - start_time))
