# Custom imports
from params        import *
from lattice_utils import *

###############################################
### LBM solver
###############################################

# Generate lattice from shape closed curve
lattice.generate_image()
lattice.solve(it_max, u_lbm, rho_lbm,
              R_ref,  U_ref, L_ref,
              output_freq)
