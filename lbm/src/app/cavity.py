# Generic imports
import math

# Custom imports
from lbm.src.app.base_app import *
from lbm.src.core.lattice import *
from lbm.src.plot.plot    import *

###############################################
### Lid-driven cavity
class cavity(base_app):
    def __init__(self):

        # Free arguments
        self.name        = 'cavity'
        self.Re_lbm      = 100.0
        self.L_lbm       = 100
        self.u_lbm       = 0.2
        self.rho_lbm     = 1.0
        self.t_max       = 20.0
        self.x_min       = 0.0
        self.x_max       = 1.0
        self.y_min       = 0.0
        self.y_max       = 1.0
        self.stop        = 'it'

        # Output parameters
        self.output_freq = 500
        self.output_it   = 0
        self.dpi         = 200

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):

        self.Cs          = 1.0/math.sqrt(3.0)
        self.ny          = self.L_lbm
        self.nu_lbm      = self.u_lbm*self.L_lbm/self.Re_lbm
        self.tau_lbm     = 0.5 + self.nu_lbm/(self.Cs**2)
        self.dt          = self.Re_lbm*self.nu_lbm/self.L_lbm**2
        self.dx          = (self.y_max-self.y_min)/self.ny
        self.dy          = self.dx
        self.nx          = math.floor(self.ny*(self.x_max-self.x_min)/
                                      (self.y_max-self.y_min))
        self.it_max      = math.floor(self.t_max/self.dt)
        self.sigma       = math.floor(10*self.nx)

    ### Initialize driven cavity fields
    def initialize(self, lattice):

        # Set fields
        self.set_inlets(lattice, 0)
        lattice.rho *= self.rho_lbm

        # Output image
        lattice.generate_image([])

        # Compute first equilibrium
        lattice.equilibrium()
        lattice.g = lattice.g_eq.copy()

    ### Set inlet fields
    def set_inlets(self, lattice, it):

        lx = lattice.lx
        ly = lattice.ly

        val  = it
        ret  = (1.0 - math.exp(-val**2/(2.0*self.sigma**2)))

        lattice.u_top[0,:]   = self.u_lbm*ret
        lattice.u_bot[0,:]   = 0.0
        lattice.u_left[1,:]  = 0.0
        lattice.u_right[1,:] = 0.0

    ### Set boundary conditions
    def set_bc(self, lattice):

        # Wall BCs
        lattice.zou_he_bottom_wall_velocity()
        lattice.zou_he_left_wall_velocity()
        lattice.zou_he_right_wall_velocity()
        lattice.zou_he_top_wall_velocity()
        lattice.zou_he_bottom_left_corner()
        lattice.zou_he_top_left_corner()
        lattice.zou_he_top_right_corner()
        lattice.zou_he_bottom_right_corner()

    ### Write outputs
    def outputs(self, lattice, it):

        # Check iteration
        if (it%self.output_freq != 0): return

        # Output field
        plot_norm(lattice, 0.0, 1.0, self.output_it, self.dpi)

        # Increment plotting counter
        self.output_it += 1

    ### Finalize
    def finalize(self, lattice):

        # Compute 1D fields to compare with ref. data
        self.line_fields(lattice)

    ### Compute line fields in the middle of the domain
    def line_fields(self, lattice):

        vx_error = np.zeros((self.nx))
        uy_error = np.zeros((self.ny))
        nx       = math.floor(self.nx/2)
        ny       = math.floor(self.ny/2)

        for i in range(self.nx):
            vx_error[i] = lattice.u[1,i,ny]/self.u_lbm

        for j in range(self.ny):
            uy_error[j] = lattice.u[0,nx,j]/self.u_lbm

        # Write to files
        filename = lattice.output_dir+'cavity_vx'
        with open(filename, 'w') as f:
            for i in range(self.nx):
                f.write('{} {}\n'.format(i*self.dx, vx_error[i]))
        filename = lattice.output_dir+'cavity_uy'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {}\n'.format(j*self.dx, uy_error[j]))

        return vx_error, uy_error
