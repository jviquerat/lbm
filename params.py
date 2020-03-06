# Generic imports
import math
import datetime

# Custom imports
from lattice_utils import *

###############################################
### LBM solver
###############################################

### Domain
x_min      = 0.0
x_max      = 1.0
y_min      = 0.0
y_max      = 1.0

### Fixed parameters
### u_lbm < 0.05 maintains low Mach condition
u_lbm      = 0.1
tau_lbm    = 0.8
dx_lbm     = 1.0
dt_lbm     = 1.0
Cs         = 1.0/math.sqrt(3.0)

### Free parameters
Re         = 400.0
rho        = 1.0
L          = x_max - x_min
nx         = 600
nu         = 1.0e-3

#nu_lbm     = u_lbm*nx/Re
nu_lbm     = 0.0001
u_lbm       = Re*nu_lbm/nx
#nx         = math.floor(Re*nu_lbm/u_lbm)

### Deduce u, dx and dt
u          = nu*Re/L
ny         = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx         = (x_max-x_min)/nx
dt         = dx*(u_lbm/u)

### Set parameters conversions
### Conversions are assumed s.t. a = Ca * a_lbm
Cx = dx
Ct = dt
Cu = Cx/Ct
Cn = Cx**2/Ct
Cr = 1.0
Cf = Cr*Cx**2/Ct

### Deduce remaining parameters
L_lbm      = L/Cx
#nu_lbm     = nu/Cn
rho_lbm    = rho/Cr
Re_lbm     = u_lbm*L_lbm/nu_lbm
tau_lbm    = 0.5 + nu_lbm/(Cs**2)
#dt         = ((tau_lbm - 0.5)*Cs**2*dx**2)/nu

### TRT parameters
tau_p_lbm  = tau_lbm
lambda_trt = 1.0/4.0 # Constant TRT parameter
tau_m_lbm  = lambda_trt/(tau_p_lbm - 0.5) + 0.5

# Output parameters
output_freq    = 500
time           = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir    = './results/'
output_dir     = results_dir+str(time)+'/'

# Other parameters
lattice_name = 'lattice'
t_max        = 30.0
it_max       = math.floor(t_max/dt) + 1
dpi          = 200

# Printings
print('### LBM solver ###')
print('# u          = '+str(u))
print('# u_lbm      = '+str(u_lbm))
print('# tau_p_lbm  = '+str(tau_p_lbm))
print('# tau_m_lbm  = '+str(tau_m_lbm))
print('# Re         = '+str(Re))
print('# Re_lbm     = '+str(Re_lbm))
print('# nx         = '+str(nx))
print('# ny         = '+str(ny))
print('# dx         = '+str(dx))
print('# dt         = '+str(dt))
print('# dx/dt      = '+str(dx/dt))
print('# nu         = '+str(nu))
print('# nu_lbm     = '+str(nu_lbm))
print('# it         = '+str(it_max))

if (not os.path.exists(results_dir)): os.makedirs(results_dir)

# Initialize lattice
lattice = Lattice(lattice_name,
                  x_min,      x_max,
                  y_min,      y_max,
                  nx,         ny,
                  tau_p_lbm,  tau_m_lbm,
                  Cx,         Ct, Cs,
                  Cr,         Cu,
                  Cf,         dx, dt,
                  output_dir, dpi)
