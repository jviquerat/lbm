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

### LBM parameters
t_max          = 1.0
Re             = 100.0
x_min          =-5.0
x_max          = 10.0
y_min          =-2.5
y_max          = 2.5
lat_density    = 100
lattice_name   = 'lattice'
nx             = math.floor((xmax-xmin)*lat_density)
ny             = math.floor((ymax-ymin)*lat_density)
it_max         = 2.0*nx

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
                  nx,    ny)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
lattice.init_computation()
