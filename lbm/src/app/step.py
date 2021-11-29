# Generic imports
import math

# Custom imports
from lbm.src.app.base_app  import *
from lbm.src.core.lattice  import *
from lbm.src.core.obstacle import *
from lbm.src.utils.buff    import *
from lbm.src.plot.plot     import *

###############################################
### A step with modified inlets/outlets
class step(base_app):
    def __init__(self):

        # Free arguments
        self.name        = 'step'
        self.Re_lbm      = 500.0
        self.L_lbm       = 150
        self.u_lbm       = 0.05
        self.rho_lbm     = 1.0
        self.t_max       = 15.0
        self.x_min       =-1.0
        self.x_max       = 15.0
        self.y_min       =-1.0
        self.y_max       = 1.0
        self.IBB         = False
        self.stop        = 'it'
        self.obs_cv_ct   = 1.0e-3
        self.obs_cv_nb   = 1000

        # Output parameters
        self.output_freq = 500
        self.output_it   = 0
        self.dpi         = 200

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

        # Obstacles
        self.obstacles = []
        square1 = obstacle('square1', 4, 100, 'square', 1.5, [1.5,1.0])
        self.obstacles.append(square1)
        square2 = obstacle('square2', 4, 100, 'square', 1.5, [4.0,-1.0])
        self.obstacles.append(square2)

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):

        self.Cs      = 1.0/math.sqrt(3.0)
        self.ny      = self.L_lbm
        self.u_avg   = 2.0*self.u_lbm/3.0
        self.r_cyl   = 0.5
        self.D_lbm   = math.floor(self.ny*self.r_cyl/(self.y_max-self.y_min))
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

        # Add obstacles to lattice
        self.add_obstacles(lattice, self.obstacles)

        # Initialize fields
        self.set_inlets(lattice, 0)
        lattice.u[:,np.where(lattice.lattice > 0.0)] = 0.0
        lattice.rho *= self.rho_lbm

        # Output image
        lattice.generate_image(self.obstacles)

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
            pt                  = lattice.get_coords(0, j)
            lattice.u_left[:,j] = ret*self.u_lbm*self.poiseuille(pt)

        lattice.u_top[0,:]   = 0.0
        lattice.u_bot[0,:]   = 0.0
        lattice.u_right[1,:] = 0.0
        lattice.rho_right[:] = self.rho_lbm

    ### Set boundary conditions
    def set_bc(self, lattice):

        # Obstacle
        lattice.bounce_back_obstacle(self.obstacles[0])
        lattice.bounce_back_obstacle(self.obstacles[1])

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

    ### Poiseuille flow
    def poiseuille(self, pt):

        x    = pt[0]
        y    = pt[1]
        H    = self.y_max - self.y_min
        u    = np.zeros(2)
        u[0] = 4.0*(self.y_max-y)*(y-self.y_min)/H**2

        return u
