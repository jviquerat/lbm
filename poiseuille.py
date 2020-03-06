# Generic imports
import progress.bar

# Custom imports
from params        import *
from lattice_utils import *

###############################################
### LBM poiseuille
###############################################

# Initialize fields
lattice.inlet_poiseuille(u_lbm)
lattice.rho *= rho_lbm

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq

# Solve
bar = progress.bar.Bar('Solving...', max=it_max)
for it in range(it_max+1):

    # Compute macroscopic fields
    lattice.macro()

    # Compute equilibrium state
    lattice.equilibrium()

    # Output view
    lattice.output_view(it, output_freq, u_lbm)

    # Collisions
    lattice.trt_collisions()

    # Streaming
    lattice.stream()

    # Boundary conditions
    lattice.zou_he_top_wall()
    lattice.zou_he_bottom_wall()
    lattice.zou_he_inlet()
    lattice.zou_he_outlet(rho_lbm)
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
