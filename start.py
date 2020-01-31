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
n_pts          = 6
n_sampling_pts = 50
shape_type     = 'cylinder' # 'cylinder' or 'random'

### LBM parameters
u_in           = 0.3
nu             = 0.001
q              = 9
x_min          =-0.2
x_max          = 2.0
y_min          =-0.2
y_max          = 0.21
U_ref          = 2.0*u_in/3.0
L_ref          = 0.1
Re             = U_ref*L_ref/nu
nx             = 500
ny             = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx             = (x_max-x_min)/nx
cs             = 1.0/math.sqrt(3.0)
tau            = 0.6
nu_lbm         = (tau - 0.5)*cs**2
dt             = (nu_lbm/nu)*dx**2
u_lbm          = u_in*(dt/dx)
U_ref         *= dt/dx
L_ref         /= dx
rho            = 1.0
lattice_name   = 'lattice'
it_max         = 50*nx
dpi            = 200

# Printings
print('### LBM solver ###')
print('# u        = '+str(u_in))
print('# u_lbm    = '+str(u_lbm)+' (should be < 0.05)')
print('# tau      = '+str(tau))
print('# Re       = '+str(Re))
print('# nx       = '+str(nx))
print('# ny       = '+str(ny))
print('# dx       = '+str(dx))
print('# dt       = '+str(dt))
print('# dx/dt    = '+str(dx/dt))
print('# it       = '+str(it_max))

# Output parameters
output_freq    = 500
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
                  q,          tau,
                  output_dir, dpi)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
lattice.solve(it_max, u_lbm, rho,
              U_ref, L_ref, dx, dt, output_freq)
