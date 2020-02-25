# Generic imports
import math
import datetime

# Custom imports
from shapes_utils  import *
from lattice_utils import *

###############################################
### LBM solver
###############################################

### Shape parameters
shape_name     = 'shape'
n_pts          = 4
n_sampling_pts = 50
shape_type     = 'cylinder' # 'cylinder', 'square' or 'random'

### Physical parameters
u_in       = 0.3
nu         = 0.001
rho        = 1.0
q          = 9
x_min      =-0.2
x_max      = 1.0
y_min      =-0.2
y_max      = 0.21
U_ref      = 2.0*u_in/3.0
L_ref      = 0.1
R_ref      = rho
Re         = U_ref*L_ref/nu
nx         = 300
ny         = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx         = (x_max-x_min)/nx
u_lim      = 0.03

### Parameters conversions
###
### All parameters in LBM units are postfixed with _lbm
### All conversion parameters   are prefixed  with C
### Conversions are assumed s.t. a = Ca * a_lbm
###
### Furthermore, we choose tau_p_lbm explicitely and deduce dt from it
Cs         = 1.0/math.sqrt(3.0)
tau_p_lbm  = 0.5 + (u_lim*nu)/(u_in*dx*Cs**2) # TRT relaxation parameters
lambda_trt = 0.25 # Constant TRT parameter
tau_m_lbm  = lambda_trt/(tau_p_lbm - 0.5) + 0.5
nu_lbm  = (tau_p_lbm - 0.5)*Cs**2
dt      = (nu_lbm/nu)*dx**2
Cx      = dx
Ct      = dt
Cr      = rho
rho_lbm = rho/Cr
Cu      = Cx/Ct
u_lbm   = u_in/Cu
#Cf      = Cr*Cx**4/Ct**2
Cf      = Cr*Cx**2/Ct

# Other parameters
lattice_name = 'lattice'
t_max        = 5.0
it_max       = math.floor(t_max/dt) + 1
dpi          = 200

# Printings
print('### LBM solver ###')
print('# u          = '+str(u_in))
print('# u_lbm      = '+str(u_lbm)+' (should be < 0.05)')
print('# tau_p_lbm  = '+str(tau_p_lbm))
print('# tau_m_lbm  = '+str(tau_m_lbm))
print('# Re         = '+str(Re))
print('# Re_lbm     = '+str(Re))
print('# nx         = '+str(nx))
print('# ny         = '+str(ny))
print('# dx         = '+str(dx))
print('# dt         = '+str(dt))
print('# dx/dt      = '+str(dx/dt))
print('# it         = '+str(it_max))

# Output parameters
output_freq    = 10
time           = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir    = './results/'
output_dir     = results_dir+str(time)+'/'

if (not os.path.exists(results_dir)): os.makedirs(results_dir)

### Generate shape
if (shape_type == 'cylinder'):
    radius         = (1.0/math.sqrt(2))*np.ones((n_pts))
    edgy           = 1.0*np.ones((n_pts))
    ctrl_pts       = generate_cylinder_pts(n_pts)
    ctrl_pts[:,:] *= 0.05

if (shape_type == 'square'):
    radius         = np.zeros((n_pts))
    edgy           = np.ones((n_pts))
    ctrl_pts       = generate_square_pts(n_pts)
    ctrl_pts[:,:] *= 0.05

if (shape_type == 'random'):
    radius         = np.random.uniform(low=0.0, high=1.0, size=n_pts)
    edgy           = np.random.uniform(low=0.0, high=1.0, size=n_pts)
    ctrl_pts       = np.random.rand(n_pts,2)

shape = Shape(shape_name,
              ctrl_pts,
              n_pts,
              n_sampling_pts,
              radius,
              edgy,
              output_dir)

shape.generate()
shape.generate_image()
shape.write_csv()

# Initialize lattice
lattice = Lattice(lattice_name,
                  x_min,      x_max,
                  y_min,      y_max,
                  nx,         ny,
                  q,          tau_p_lbm, tau_m_lbm,
                  Cx,         Ct,
                  Cr,         Cu,
                  Cf,
                  output_dir, dpi)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
lattice.solve(it_max, u_lbm, rho_lbm,
              R_ref,  U_ref, L_ref,
              output_freq)
