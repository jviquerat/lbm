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
q              = 9     # D2Q9 lattice
x_min          =-5.0
x_max          = 10.0
y_min          =-2.5
y_max          = 2.5
lat_density    = 50
lattice_name   = 'lattice'
nx             = math.floor((x_max-x_min)*lat_density)
ny             = math.floor((y_max-y_min)*lat_density)
it_max         = 5*nx

### Fluid parameters
Re             = 100.0
u_in           = 0.04
nu             = u_in*ny/Re
tau            = 0.5 + 3.0*nu # relaxation parameter
rho            = 1.0       # fluid density

# Output parameters
output_freq    = 100
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
lattice.init_computation(u_in)
lattice.solve(it_max, u_in, rho, output_freq)
