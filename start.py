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
radius         = np.random.uniform(low=0.0, high=1.0, size=n_pts)
edgy           = np.random.uniform(low=0.0, high=1.0, size=n_pts)
ctrl_pts       = np.random.rand(n_pts,2)

### LBM parameters
u              = 0.1
nu             = 0.1
q              = 9
x_min          =-5.0
x_max          = 10.0
y_min          =-5.0
y_max          = 5.0
Re             = u*(x_max-x_min)/nu
nx             = 300
ny             = math.floor(nx*(y_max-y_min)/(x_max-x_min))
dx             = (x_max-x_min)/nx
cs             = 1.0/math.sqrt(3.0)
tau            = 0.55
nu_lbm         = (tau - 0.5)*cs**2
dt             = (nu_lbm/nu)*dx**2
u_lbm          = u/(dx/dt)
rho            = 1.0

#lat_density    = 25
lattice_name   = 'lattice'
#nx             = math.floor((x_max-x_min)*lat_density)
#ny             = math.floor((y_max-y_min)*lat_density)
it_max         = 20*nx

### Fluid parameters
#Re             = 100.0  # Reynolds number
#cs             = 1.0/math.sqrt(3.0)  # Speed of sound
#L              = ny # typical size
#dx             = L/ny
#dt             = dx/cs
#u_in           = dx/dt
#nu             = u_in*L/Re
#tau            = 0.5 + 3.0*nu # relaxation parameter


print('### LBM solver ###')
print('# u        = '+str(u))
print('# u_lbm    = '+str(u_lbm))
print('# tau      = '+str(tau))
print('# u_lbm/cs = '+str(u_lbm/cs))
print('# Re       = '+str(Re))
print('# nx       = '+str(nx))
print('# ny       = '+str(ny))
print('# dx       = '+str(dx))
print('# dt       = '+str(dt))
print('# it       = '+str(it_max))

# Output parameters
output_freq    = 10
time           = datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
results_dir    = './results/'
output_dir     = results_dir+str(time)+'/'

if (not os.path.exists(results_dir)): os.makedirs(results_dir)

### Initialize shape
shape = Shape(shape_name,
              ctrl_pts,
              n_pts,
              n_sampling_pts,
              radius,
              edgy)

# Generate shape
shape.generate()
shape.generate_image()
shape.write_csv()

# Initialize lattice
lattice = Lattice(lattice_name,
                  x_min, x_max,
                  y_min, y_max,
                  nx,    ny,
                  q,     tau,
                  output_dir)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
lattice.init_computation()
lattice.solve(it_max, u_lbm, rho, output_freq)
