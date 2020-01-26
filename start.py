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
u              = 1.0
nu             = 0.01
q              = 9
x_min          =-5.0
x_max          = 15.0
y_min          =-2.0
y_max          = 2.0
Re             = u*(y_max-y_min)/nu
nx             = 1000
ny             = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx             = (x_max-x_min)/nx
cs             = 1.0/math.sqrt(3.0)
tau            = 0.55
nu_lbm         = (tau - 0.5)*cs**2
dt             = (nu_lbm/nu)*dx**2
u_lbm          = u*(dt/dx)
rho            = 1.0
lattice_name   = 'lattice'
it_max         = 20*nx
dpi            = 200

# Printings
print('### LBM solver ###')
print('# u        = '+str(u))
print('# u_lbm    = '+str(u_lbm)+' (should be < 0.05)')
print('# tau      = '+str(tau))
print('# Re       = '+str(Re))
print('# nx       = '+str(nx))
print('# ny       = '+str(ny))
print('# dx       = '+str(dx))
print('# dt       = '+str(dt))
print('# it       = '+str(it_max))

# Output parameters
output_freq    = 50
time           = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir    = './results/'
output_dir     = results_dir+str(time)+'/'

if (not os.path.exists(results_dir)): os.makedirs(results_dir)

### Generate shape
if (shape_type == 'cylinder'):
    radius         = (1.0/math.sqrt(2))*np.ones((n_pts))
    edgy           = 1.0*np.ones((n_pts))
    ctrl_pts       = generate_cylinder_pts(n_pts)
    ctrl_pts[:,1] += 0.25

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
lattice.init_computation()
lattice.solve(it_max, u_lbm, rho, output_freq)
