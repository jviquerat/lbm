# Generic imports
import math

###############################################
### Base app
class base_app():
    def __init__(self):

        # Default free arguments (definition is mandatory)
        self.Re_lbm      = 100.0
        self.L_lbm       = 200
        self.u_lbm       = 0.03
        self.rho_lbm     = 1.0
        self.t_max       = 20.0

        # Output parameters
        self.output_freq = 500
        self.dpi         = 200

    ### Compute default lbm parameters
    def compute_lbm_parameters(self):

        self.Cs          = 1.0/math.sqrt(3.0)
        self.nx          = self.L_lbm
        self.ny          = self.L_lbm
        self.nu_lbm      = self.u_lbm*self.L_lbm/self.Re_lbm
        self.tau_lbm     = 0.5 + self.nu_lbm/(self.Cs**2)
        self.dt          = self.Re_lbm*self.nu_lbm/self.L_lbm**2
        self.dx          = 1.0/float(self.nx)
        self.dy          = 1.0/float(self.ny)
        self.it_max      = math.floor(self.t_max/self.dt)
