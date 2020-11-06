# Generic imports
import math
import time

# Custom imports
from lattice import *

###############################################
### LBM poiseuille
###############################################

# Domain size
x_min       =-0.2
x_max       = 2.0
y_min       =-0.2
y_max       = 0.2

# Free parameters
# L_lbm corresponds to y   length
# u_lbm corresponds to max velocity
Re_lbm      = 100.0
u_lbm       = 0.1
L_lbm       = 100
t_max       = 8.0

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
ny          = L_lbm
u_avg       = u_lbm/2.0
nu_lbm      = u_avg*L_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
it_max      = math.floor(t_max/dt)
dx          = (y_max-y_min)/ny
dy          = dx
nx          = math.floor(ny*(x_max-x_min)/(y_max-y_min))

# Other parameters
output_freq = 500
dpi         = 200
sigma       = math.floor(it_max/5)

# Initialize lattice
lattice = Lattice(nx      = nx,
                  ny      = ny,
                  dx      = dx,
                  dt      = dt,
                  tau_lbm = tau_lbm,
                  Re_lbm  = Re_lbm,
                  u_lbm   = u_lbm,
                  L_lbm   = L_lbm,
                  nu_lbm  = nu_lbm,
                  rho_lbm = rho_lbm,
                  x_min   = x_min,
                  x_max   = x_max,
                  y_min   = y_min,
                  y_max   = y_max,
                  dpi     = dpi,
                  t_max   = t_max,
                  it_max  = it_max)

# Initialize fields
lattice.set_inlet_poiseuille(u_lbm, rho_lbm, 0, sigma)

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq.copy()

# Count time
start_time = time.time()

# Solve
print('### Solving')
while (lattice.compute):

    # Printings
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

    # Streaming
    lattice.collision_stream()

    # Boundary conditions
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_right_wall_pressure()
    lattice.zou_he_top_wall_velocity()
    lattice.zou_he_bottom_left_corner()
    lattice.zou_he_top_left_corner()
    lattice.zou_he_top_right_corner()
    lattice.zou_he_bottom_right_corner()

    # Check stopping criterion
    lattice.check_stop()

# Count time
end_time = time.time()
print("# Loop time = {:f}".format(end_time - start_time))

# Output error with exact solution
lattice.poiseuille_error(u_lbm)
