# Generic imports
import math
import progress.bar

# Custom imports
from lattice_utils import *

###############################################
### LBM lid-driven cavity
###############################################

### Free parameters
Re_lbm      = 100.0
u_lbm       = 0.03
L_lbm       = 200#600

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
nx          = L_lbm
nu_lbm      = u_lbm*L_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0

# Other parameters
output_freq = 500
dpi         = 200
it_max      = 120000

# Initialize lattice
lattice = Lattice(nx      = nx,
                  dx      = 1.0/nx,
                  tau_lbm = tau_lbm,
                  Re_lbm  = Re_lbm,
                  u_lbm   = u_lbm,
                  L_lbm   = L_lbm,
                  nu_lbm  = nu_lbm,
                  rho_lbm = rho_lbm,
                  dpi     = dpi)

# Initialize fields
lattice.set_cavity(u_lbm)
lattice.rho *= rho_lbm

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq.copy()

# Solve
bar = progress.bar.Bar('Solving...', max=it_max)
for it in range(it_max+1):

    # Compute macroscopic fields
    lattice.macro()

    # Boundary conditions
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_right_wall_velocity()
    lattice.zou_he_top_wall_velocity()
    lattice.zou_he_bottom_left_corner()
    lattice.zou_he_top_left_corner()
    lattice.zou_he_top_right_corner()
    lattice.zou_he_bottom_right_corner()

    # Compute equilibrium state
    lattice.equilibrium()

    # Collisions
    lattice.trt_collisions()

    # Streaming
    lattice.stream()

    # Output view
    lattice.output_fields(it,
                          output_freq,
                          u_norm   = True,
                          u_stream = True)

    # Increment bar
    bar.next()

# End bar
bar.finish()

# Output error with exact solution
lattice.cavity_error(u_lbm)
