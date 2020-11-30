# Generic imports
import math
import time

# Custom imports
from lattice import *

###############################################
### LBM lid-driven cavity
###############################################

### Free parameters
Re_lbm      = 100.0
u_lbm       = 0.03
L_lbm       = 400
t_max       = 8.0

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
nx          = L_lbm
nu_lbm      = u_lbm*L_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
it_max      = math.floor(t_max/dt)

# Other parameters
output_freq = 2000
dpi         = 200

# Initialize lattice
lattice = Lattice(nx      = nx,
                  ny      = nx,
                  dx      = 1.0/nx,
                  dt      = dt,
                  tau_lbm = tau_lbm,
                  Re_lbm  = Re_lbm,
                  u_lbm   = u_lbm,
                  L_lbm   = L_lbm,
                  nu_lbm  = nu_lbm,
                  rho_lbm = rho_lbm,
                  dpi     = dpi,
                  t_max   = t_max,
                  it_max  = it_max)

# Initialize fields
lattice.set_cavity(-u_lbm, -2.0*u_lbm, -1.0*u_lbm, 1.5*u_lbm)
lattice.rho *= rho_lbm

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

    # Compute macroscopic fields
    lattice.macro()

    # Output field
    lattice.output_fields(lattice.it,
                          output_freq,
                          u_norm   = True,
                          u_ctr    = False,
                          u_stream = False)

    # Compute equilibrium state
    lattice.equilibrium()

    # Streaming
    lattice.collision_stream()

    # Boundary conditions
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_right_wall_velocity()
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
lattice.cavity_error(u_lbm)

# Output streamlines
lattice.output_fields(1, 1,
                      u_norm   = False,
                      u_stream = True)
