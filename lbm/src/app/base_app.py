# Generic imports
import math

# Custom imports
from lbm.src.utils.shapes import *
from lbm.src.core.obstacle import *

###############################################
### Base app
class base_app():
    def __init__(self):

        # Default free arguments (definition is mandatory)
        self.name        = 'base'
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.03
        self.rho_lbm     = 1.0
        self.t_max       = 20.0
        self.x_min       = 0.0
        self.x_max       = 1.0
        self.y_min       = 0.0
        self.y_max       = 1.0

        # Output parameters
        self.output_freq = 500
        self.output_it   = 0
        self.dpi         = 200

    ### Compute default lbm parameters
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

    ### Set inlet fields
    def set_inlets(self, lattice, it):

        pass

    ### Compute observables
    def observables(self, lattice, it):

        pass

    ### Add obstacles
    def add_obstacles(self, lattice, obstacles):

        for i in range(len(obstacles)):
            obs   = obstacles[i]
            shape = generate_shape(obs.n_pts, obs.pos,
                                   obs.type,  obs.size,
                                   obs.name,  obs.n_spts,
                                   lattice.output_dir)
            obstacles[i].set_polygon(shape.curve_pts)
            obstacles[i].set_tag(i)
            area, bnd, ibb = lattice.add_obstacle(obstacles[i])
            obstacles[i].fill(area, bnd, ibb)
