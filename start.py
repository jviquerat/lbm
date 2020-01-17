# Generic imports
import math

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

### Fluid parameters
u_in           = 1.0       # input velocity
rho            = 1.0       # fluid density
mu             = 0.1       # dynamic viscosity
nu             = mu/rho    # cinematic viscosity
L              = 1.0       # baseline dimension
Re             = u_in*L/nu # baseline Reynolds for L=1

### LBM parameters
q              = 9     # D2Q9 lattice
t_max          = 1.0

x_min          =-5.0
x_max          = 10.0
y_min          =-2.5
y_max          = 2.5
lat_density    = 10
lattice_name   = 'lattice'
nx             = math.floor((x_max-x_min)*lat_density)
ny             = math.floor((y_max-y_min)*lat_density)
it_max         = 2*nx

### Normalize units
c_sound        = 343.0 # speed of sound
u_in           = Re/c_sound


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
                  q)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
lattice.init_computation(u_in)
lattice.solve(it_max)
