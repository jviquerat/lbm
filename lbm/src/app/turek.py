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

        # Free arguments
        self.name        = 'turek'
        self.Re_lbm      = 100.0
        self.L_lbm       = 151
        self.u_lbm       = 0.05
        self.rho_lbm     = 1.0
        self.t_max       = 0.02
        self.x_min       =-0.2
        self.x_max       = 2.0
        self.y_min       =-0.2
        self.y_max       = 0.21
        self.IBB         = True
        self.stop        = 'obs'
        self.obs_cv_ct   = 1.0e-3
        self.obs_cv_nb   = 1000

        # Output parameters
        self.output_freq = 500
        self.output_it   = 0
        self.dpi         = 200

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

        # Obstacle
        cylinder = obstacle('turek', 200, 2,
                            'cylinder', self.r_cyl, [0.0,0.0])
        self.obstacles = [cylinder]

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):

        self.Cs      = 1.0/math.sqrt(3.0)
        self.ny      = self.L_lbm
        self.u_avg   = 2.0*self.u_lbm/3.0
        self.r_cyl   = 0.1
        self.D_lbm   = math.floor(self.ny*self.r_cyl/(self.y_max-self.y_min))
        self.nu_lbm  = self.u_avg*self.D_lbm/self.Re_lbm
        self.tau_lbm = 0.5 + self.nu_lbm/(self.Cs**2)
        self.dt      = self.Re_lbm*self.nu_lbm/self.D_lbm**2
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

        # Obstacle
        lattice.bounce_back_obstacle(self.obstacles[0])

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

    ### Compute observables
    def observables(self, lattice, it):

        drag, lift = lattice.drag_lift(self.obstacles[0],
                                       self.rho_lbm, self.u_avg, self.D_lbm)
        self.add_buff(drag, lift, lattice, it)

    ### Handle drag/lift buffers
    def add_buff(self, Cx, Cy, lattice, it):

        # Add to buffer and check for convergence
        self.drag_buff.add(Cx)
        self.lift_buff.add(Cy)

        avg_Cx, dcx = self.drag_buff.mv_avg()
        avg_Cy, dcy = self.lift_buff.mv_avg()

        # Write to file
        filename = lattice.output_dir+'drag_lift'
        with open(filename, 'a') as f:
            f.write('{} {} {} {} {} {} {}\n'.format(it*self.dt,
                                                    Cx,     Cy,
                                                    avg_Cx, avg_Cy,
                                                    dcx,    dcy))

    ### Poiseuille flow
    def poiseuille(self, pt):

        x    = pt[0]
        y    = pt[1]
        H    = self.y_max - self.y_min
        u    = np.zeros(2)
        u[0] = 4.0*(self.y_max-y)*(y-self.y_min)/H**2

        return u
