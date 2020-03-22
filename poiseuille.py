# Generic imports
import math
import progress.bar

# Custom imports
from lattice_utils import *

###############################################
### LBM poiseuille
###############################################

# Domain size
x_min       =-0.2
x_max       = 1.0
y_min       =-0.2
y_max       = 0.2

# Free parameters
# L_lbm corresponds to y   length
# u_lbm corresponds to max velocity
Re_lbm      = 100.0
u_lbm       = 0.1
L_lbm       = 100
t_max       = 15.0

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

# Poiseuille imposition time
sigma       = math.floor(it_max/10)

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
                  dpi     = dpi)

# Initialize fields
lattice.set_poiseuille(u_lbm, rho_lbm, 0, sigma)

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq.copy()

# Solve
bar = progress.bar.Bar('Solving...', max=it_max)
for it in range(it_max+1):

    lattice.set_poiseuille(u_lbm, rho_lbm, it, sigma)

    # Compute macroscopic fields
    lattice.macro()

    # Output field
    lattice.output_fields(it,
                          output_freq,
                          u_norm   = True,
                          u_stream = False)

    # Compute equilibrium state
    lattice.equilibrium()

    # Collisions
    lattice.trt_collisions()

    # Streaming
    lattice.stream()

    # Compute equilibrium state
    lattice.equilibrium()

    # Boundary conditions
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_right_wall_pressure()
    lattice.zou_he_top_wall_velocity()
    lattice.zou_he_bottom_left_corner()
    lattice.zou_he_top_left_corner()
    lattice.zou_he_top_right_corner()
    lattice.zou_he_bottom_right_corner()

    # Increment bar
    bar.next()

# End bar
bar.finish()

# Output error with exact solution
lattice.poiseuille_error(u_lbm)
