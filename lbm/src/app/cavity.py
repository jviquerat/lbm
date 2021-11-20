# Generic imports
import math

# Custom imports
from lbm.src.app.base_app import *
from lbm.src.core.lattice import *
from lbm.src.plot.plot    import *

###############################################
### Lid-drive cavity
class cavity(base_app):
    def __init__(self):
        super().__init__()

        # Free arguments (definition is mandatory)
        self.name        = 'cavity'
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.1
        self.rho_lbm     = 1.0
        self.t_max       = 12.0
        self.x_min       = 0.0
        self.x_max       = 1.0
        self.y_min       = 0.0
        self.y_max       = 1.0

        # Output parameters
        self.output_freq = 5000

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

    ### Initialize driven cavity fields
    def initialize(self, lattice):

        self.top_bc(lattice, 0)
        lattice.rho *= self.rho_lbm

    ### Set top bc
    def top_bc(self, lattice, it):

        lx = lattice.lx
        ly = lattice.ly

        lattice.u_top[0,:]   = self.u_lbm
        lattice.u_bot[0,:]   = 0.0
        lattice.u_left[1,:]  = 0.0
        lattice.u_right[1,:] = 0.0

        lattice.u[:,:,ly]  = lattice.u_top[:,:]
        lattice.u[:,:,0]   = 0.0
        lattice.u[:,0,:]   = 0.0
        lattice.u[:,lx,:]  = 0.0

    ### Write outputs
    def outputs(self, lattice, it):

        # Check iteration
        if (it%self.output_freq != 0): return

        # Output field
        plot_norm(lattice, self.output_it, self.dpi)

        # Increment plotting counter
        self.output_it += 1

    ### Finalize
    def finalize(self, lattice):

        # Compute 1D fields to compare with ref. data
        self.line_fields(lattice)

        # Write streamline field
        #plot_streamlines(lattice, self.dpi)

    ### Compute line fields in the middle of the domain
    def line_fields(self, lattice):

        ux_error = np.zeros((self.nx))
        uy_error = np.zeros((self.ny))
        nx       = math.floor(self.nx/2)
        ny       = math.floor(self.ny/2)

        for i in range(self.nx):
            uy_error[i] = lattice.u[1,i,ny]/self.u_lbm

        for j in range(self.ny):
            ux_error[j] = lattice.u[0,nx,j]/self.u_lbm

        # Write to files
        filename = lattice.output_dir+'cavity_uy'
        with open(filename, 'w') as f:
            for i in range(self.nx):
                f.write('{} {}\n'.format(i*self.dx, uy_error[i]))
        filename = lattice.output_dir+'cavity_ux'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {}\n'.format(j*self.dx, ux_error[j]))
