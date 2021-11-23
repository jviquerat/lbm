# Generic imports
import math

# Custom imports
from lbm.src.app.base_app  import *
from lbm.src.core.lattice  import *
from lbm.src.core.obstacle import *
from lbm.src.utils.buff    import *
from lbm.src.plot.plot     import *

###############################################
### Turek benchmark
class turek(base_app):
    def __init__(self):
        super().__init__()

        # Free arguments (definition is mandatory)
        self.name        = 'turek'
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.2
        self.rho_lbm     = 1.0
        self.t_max       = 0.02
        self.x_min       =-0.2
        self.x_max       = 2.0
        self.y_min       =-0.2
        self.y_max       = 0.21
        self.IBB         = True
        self.stop        = 'obs'
        self.obs_cv_ct   = 1.0e-2
        self.obs_cv_nb   = 1000

        # Output parameters
        self.output_freq = 10000

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

        # Obstacle
        cylinder = obstacle('turek', 200, 2, 'cylinder', 0.1, [0.0,0.0])
        self.obstacles = [cylinder]

    ### Add obstacles and initialize fields
    def initialize(self, lattice):

        # Add obstacles to lattice
        self.add_obstacles(lattice, self.obstacles)

        # Initialize fields
        self.set_bc(lattice, 0)
        lattice.u[:,np.where(lattice.lattice > 0.0)] = 0.0
        lattice.rho *= self.rho_lbm

        # Output image
        lattice.generate_image(self.obstacles)

        # Set buffers
        self.drag_buff = buff('drag',
                              lattice.dt,
                              lattice.obs_cv_ct,
                              lattice.obs_cv_nb,
                              lattice.output_dir)
        self.lift_buff = buff('lift',
                              lattice.dt,
                              lattice.obs_cv_ct,
                              lattice.obs_cv_nb,
                              lattice.output_dir)

        # Compute first equilibrium
        lattice.equilibrium()
        lattice.g = lattice.g_eq.copy()

    ### Set inlet fields
    def set_inlets(self, lattice, it):

        self.set_bc(lattice, it)

    ### Set inlet poiseuille
    def set_bc(self, lattice, it):

        lx = lattice.lx
        ly = lattice.ly

        val  = it
        ret  = (1.0 - math.exp(-val**2/(2.0*self.sigma**2)))

        for j in range(self.ny):
            pt               = lattice.lattice_coords(0, j)
            lattice.u_left[:,j] = ret*self.u_lbm*self.poiseuille(pt)

        lattice.u_top[0,:]   = 0.0
        lattice.u_bot[0,:]   = 0.0
        lattice.u_right[1,:] = 0.0
        lattice.rho_right[:] = self.rho_lbm

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

    ### Poiseuille flow
    def poiseuille(self, pt):

        x    = pt[0]
        y    = pt[1]
        H    = self.y_max - self.y_min
        u    = np.zeros(2)
        u[0] = 4.0*(self.y_max-y)*(y-self.y_min)/H**2

        return u
