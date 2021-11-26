# Generic imports
import math

# Custom imports
from lbm.src.app.base_app  import *
from lbm.src.core.lattice  import *
from lbm.src.core.obstacle import *
from lbm.src.utils.buff    import *
from lbm.src.plot.plot     import *

###############################################
### Poiseuille benchmark
class poiseuille(base_app):
    def __init__(self):

        # Free arguments
        self.name        = 'poiseuille'
        self.Re_lbm      = 100.0
        self.L_lbm       = 50
        self.u_lbm       = 0.1
        self.rho_lbm     = 1.0
        self.t_max       = 15.0
        self.x_min       =-0.2
        self.x_max       = 1.0
        self.y_min       =-0.2
        self.y_max       = 0.2
        self.stop        = 'it'

        # Output parameters
        self.output_freq = 50000
        self.output_it   = 0
        self.dpi         = 200

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):

        self.Cs      = 1.0/math.sqrt(3.0)
        self.ny      = self.L_lbm
        self.u_avg   = 2.0*self.u_lbm/3.0
        self.nu_lbm  = self.u_avg*self.L_lbm/self.Re_lbm
        self.tau_lbm = 0.5 + self.nu_lbm/(self.Cs**2)
        self.dt      = self.Re_lbm*self.nu_lbm/self.L_lbm**2
        self.dx      = (self.y_max-self.y_min)/self.ny
        self.dy      = self.dx
        self.nx      = math.floor(self.ny*(self.x_max-self.x_min)/
                                  (self.y_max-self.y_min))
        self.it_max  = math.floor(self.t_max/self.dt)
        self.sigma   = math.floor(10*self.nx)

    ### Add obstacles and initialize fields
    def initialize(self, lattice):

        # Initialize fields
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

        for j in range(self.ny):
            pt               = lattice.get_coords(0, j)
            lattice.u_left[:,j] = ret*self.u_lbm*self.poiseuille(pt)

        lattice.u_top[0,:]   = 0.0
        lattice.u_bot[0,:]   = 0.0
        lattice.u_right[1,:] = 0.0
        lattice.rho_right[:] = self.rho_lbm

    ### Set boundary conditions
    def set_bc(self, lattice):

        # Wall BCs
        lattice.zou_he_bottom_wall_velocity()
        lattice.zou_he_left_wall_velocity()
        lattice.zou_he_top_wall_velocity()
        lattice.zou_he_right_wall_pressure()
        lattice.zou_he_bottom_left_corner()
        lattice.zou_he_top_left_corner()
        lattice.zou_he_top_right_corner()
        lattice.zou_he_bottom_right_corner()

    ### Write outputs
    def outputs(self, lattice, it):

        # Check iteration
        if (it%self.output_freq != 0): return

        # Output field
        plot_norm(lattice, 0.0, 1.5, self.output_it, self.dpi)

        # Increment plotting counter
        self.output_it += 1

    ### Finalize
    def finalize(self, lattice):

        # Compute error
        self.compute_error(lattice)

    ### Poiseuille flow
    def poiseuille(self, pt):

        x    = pt[0]
        y    = pt[1]
        H    = self.y_max - self.y_min
        u    = np.zeros(2)
        u[0] = 4.0*(self.y_max-y)*(y-self.y_min)/H**2

        return u

    ### Poiseuille error in the middle of the domain
    def compute_error(self, lattice):

        u_error  = np.zeros((3,self.ny))
        nx       = math.floor(self.nx/2)
        l1_error = 0.0

        for j in range(self.ny):
            pt   = lattice.get_coords(nx,j)
            u_ex = self.poiseuille(pt)
            u    = lattice.u[:,nx,j]

            u_error[0,j] = u[0]/self.u_lbm
            u_error[1,j] = u_ex[0]
            u_error[2,j] = abs(u_error[0,j] - u_error[1,j])
            l1_error    += u_error[2,j]*self.dx

        # Write to file
        filename = lattice.output_dir+'poiseuille'
        with open(filename, 'w') as f:
            for j in range(self.ny):
                f.write('{} {} {}\n'.format(j*self.dx,
                                            u_error[0,j],
                                            u_error[1,j]))

        return l1_error
