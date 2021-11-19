# Generic imports
import math

# Custom imports
from lbm.src.app.base_app import *
from lbm.src.core.lattice import *

###############################################
### Lid-drive cavity
class cavity(base_app):
    def __init__(self):

        # Free arguments (definition is mandatory)
        self.name        = 'cavity'
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.03
        self.rho_lbm     = 1.0
        self.t_max       = 20.0
        self.x_min       = 0.0
        self.x_max       = 1.0
        self.y_min       = 0.0
        self.y_max       = 1.0

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()



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

    ### Initialize driven cavity fields
    def init_fields(self, lattice):

        self.inlet_fields(lattice, 0)
        lattice.rho *= self.rho_lbm

    ### Set inlet fields
    def inlet_fields(self, lattice, it):

        lx = lattice.lx
        ly = lattice.ly

        lattice.u_top[0,:] = self.u_lbm

        lattice.u[:,:,ly]  = lattice.u_top[:,:]
        lattice.u[:,:,0]   = 0.0
        lattice.u[:,0,:]   = 0.0
        lattice.u[:,lx,:]  = 0.0

        # self.u[0,:,ly]   = self.u_top[0,:]
        # self.u[1,:,ly]   = self.u_top[1,:]
        # self.u[0,:,0]    = 0.0
        # self.u[1,:,0]    = 0.0
        # self.u[0,0,:]    = 0.0
        # self.u[1,0,:]    = 0.0
        # self.u[0,lx,:]   = 0.0
        # self.u[1,lx,:]   = 0.0

    ### Compute cavity error in the middle of the domain
    def compute_error(self, u):

        ux_error = np.zeros((self.nx))
        uy_error = np.zeros((self.ny))
        nx       = math.floor(self.nx/2)
        ny       = math.floor(self.ny/2)

        for i in range(self.nx):
            uy_error[i] = u[1,i,ny]/self.u_lbm

        for j in range(self.ny):
            ux_error[j] = u[0,nx,j]/self.u_lbm

        # Write to files
        filename = self.output_dir+'cavity_uy'
        with open(filename, 'w') as f:
            for i in range(self.nx):
                f.write('{} {}\n'.format(i*self.dx, uy_error[i]))
        filename = self.output_dir+'cavity_ux'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {}\n'.format(j*self.dx, ux_error[j]))
