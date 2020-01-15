# Custom imports
from shapes_utils  import *
from lattice_utils import *

###############################################
### 
###############################################

### Shape parameters
shape_name     = 'shape'
n_pts          = 6
n_sampling_pts = 50
radius         = np.random.uniform(low=0.0, high=1.0, size=n_pts)
edgy           = np.random.uniform(low=0.0, high=1.0, size=n_pts)
ctrl_pts       = np.random.rand(n_pts,2)

### LBM parameters
lattice_name   = 'lattice'
Re             = 100.0
xmin           =-10.0
xmax           = 20.0
ymin           =-5.0
ymax           = 5.0
nx             = 300
ny             = 300

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

# Initialize lattice
lattice = Lattice(lattice_name,
                  xmin, xmax,
                  ymin, ymax,
                  nx,   ny)

# Generate lattice from shape closed curve
lattice.generate(shape.curve_pts)
lattice.generate_image()
