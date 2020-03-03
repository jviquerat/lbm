# Generic imports
import math
import datetime

# Custom imports
from lattice_utils import *

###############################################
### LBM solver
###############################################

### Domain
x_min      =-0.2
x_max      = 1.0
y_min      =-0.2
y_max      = 0.2

### Free parameters
nx         = 200
u_lbm      = 0.05
Re         = 10.0
nu         = 0.001
rho        = 1.0

### Deduce remaining parameters
L_ref      = 0.1
R_ref      = rho
ny         = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx         = (x_max-x_min)/nx

### Parameters conversions
###
### All parameters in LBM units are postfixed with _lbm
### All conversion parameters   are prefixed  with C
### Conversions are assumed s.t. a = Ca * a_lbm
Cx      = dx
Cr      = rho
Cs      = 1.0/math.sqrt(3.0)

L_lbm   = L_ref/Cx
U_ref   = nu*Re/L_ref
u_in    = 3.0*U_ref/2.0
nu_lbm  = u_lbm*L_lbm/Re
Cnu     = nu/nu_lbm
Ct      = Cx**2/Cnu
dt      = Ct
rho_lbm = rho/Cr
Cu      = Cx/Ct
Cf      = Cr*Cx**2/Ct
Re_lbm  = u_lbm*L_lbm/nu_lbm

tau_p_lbm  = 0.5 + nu_lbm/(Cs**2)
lambda_trt = 1.0/4.0 # Constant TRT parameter
tau_m_lbm  = lambda_trt/(tau_p_lbm - 0.5) + 0.5

# Output parameters
output_freq    = 500
time           = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir    = './results/'
output_dir     = results_dir+str(time)+'/'

# Other parameters
lattice_name = 'lattice'
t_max        = 100.0
it_max       = math.floor(t_max/dt) + 1
dpi          = 200

# Printings
print('### LBM solver ###')
print('# u          = '+str(u_in))
print('# u_lbm      = '+str(u_lbm)+' (should be < 0.05)')
print('# tau_p_lbm  = '+str(tau_p_lbm))
print('# tau_m_lbm  = '+str(tau_m_lbm))
print('# Re         = '+str(Re))
print('# Re_lbm     = '+str(Re_lbm))
print('# nx         = '+str(nx))
print('# ny         = '+str(ny))
print('# dx         = '+str(dx))
print('# dt         = '+str(dt))
print('# dx/dt      = '+str(dx/dt))
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
