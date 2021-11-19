# Generic imports
import math
import time

# Custom imports
from lbm.src.app.base_app import *
from lbm.src.core.lattice import *

###############################################
### Lid-drive cavity
class cavity(base_app):
    def __init__(self):

        # Free arguments (definition is mandatory)
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.03
        self.rho_lbm     = 1.0
        self.t_max       = 20.0

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()


# # Initialize lattice
# lattice = Lattice(nx      = nx,
#                   ny      = nx,
#                   dx      = 1.0/nx,
#                   dt      = dt,
#                   tau_lbm = tau_lbm,
#                   Re_lbm  = Re_lbm,
#                   u_lbm   = u_lbm,
#                   L_lbm   = L_lbm,
#                   nu_lbm  = nu_lbm,
#                   rho_lbm = rho_lbm,
#                   dpi     = dpi,
#                   t_max   = t_max,
#                   it_max  = it_max)

# # Initialize fields
# lattice.set_cavity(u_lbm)
# lattice.rho *= rho_lbm

# # Set initial distributions
# lattice.equilibrium()
# lattice.g = lattice.g_eq.copy()

# # Count time
# start_time = time.time()

# # Solve
# print('### Solving')
# while (lattice.compute):

#     # Printings
#     lattice.it_printings()

#     # Compute macroscopic fields
#     lattice.macro()

#     # Output field
#     lattice.output_fields(lattice.it,
#                           output_freq,
#                           u_norm   = True,
#                           u_stream = False)

#     # Compute equilibrium state
#     lattice.equilibrium()

#     # Streaming
#     lattice.collision_stream()

#     # Boundary conditions
#     lattice.zou_he_bottom_wall_velocity()
#     lattice.zou_he_left_wall_velocity()
#     lattice.zou_he_right_wall_velocity()
#     lattice.zou_he_top_wall_velocity()
#     lattice.zou_he_bottom_left_corner()
#     lattice.zou_he_top_left_corner()
#     lattice.zou_he_top_right_corner()
#     lattice.zou_he_bottom_right_corner()

#     # Check stopping criterion
#     lattice.check_stop()

# # Count time
# end_time = time.time()
# print("# Loop time = {:f}".format(end_time - start_time))

# # Output error with exact solution
# lattice.cavity_error(u_lbm)

# # Output streamlines
# lattice.output_fields(1, 1,
#                       u_norm   = False,
#                       u_stream = True)
