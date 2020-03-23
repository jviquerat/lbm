# Generic imports
import math
import progress.bar

# Custom imports
from shapes_utils  import *
from lattice_utils import *

###############################################
### LBM poiseuille with obstacle
###############################################

# Shape parameters
shape_name     = 'shape'
n_pts          = 4
n_sampling_pts = 50
shape_type     = 'cylinder' # 'cylinder', 'square' or 'random'
shape_size     = 0.1

# Domain size
x_min       =-0.2
x_max       = 2.0
y_min       =-0.2
y_max       = 0.21

# Free parameters
# L_lbm correponds to y length
# u_lbm corresponds to max velocity
Re_lbm      = 20.0
u_lbm       = 0.025
L_lbm       = 50
t_max       = 1.0

# Deduce other parameters
Cs          = 1.0/math.sqrt(3.0)
ny          = L_lbm
u_avg       = 2.0*u_lbm/3.0
D_lbm       = math.floor(ny*shape_size/(y_max-y_min))
nu_lbm      = u_avg*D_lbm/Re_lbm
tau_lbm     = 0.5 + nu_lbm/(Cs**2)
rho_lbm     = 1.0
dt          = Re_lbm*nu_lbm/L_lbm**2
it_max      = math.floor(t_max/dt)
dx          = (y_max-y_min)/ny
dy          = dx
nx          = math.floor(ny*(x_max-x_min)/(y_max-y_min))

# Other parameters
output_freq = 500
dpi         = 200

# Initialize lattice
lattice = Lattice(nx      = nx,
                  ny      = ny,
                  dx      = dx,
                  dt      = dt,
                  tau_lbm = tau_lbm,
                  Re_lbm  = Re_lbm,
                  u_lbm   = u_lbm,
                  L_lbm   = D_lbm,
                  nu_lbm  = nu_lbm,
                  rho_lbm = rho_lbm,
                  x_min   = x_min,
                  x_max   = x_max,
                  y_min   = y_min,
                  y_max   = y_max,
                  dpi     = dpi)

# Generate shape
if (shape_type == 'cylinder'):
    radius         = (1.0/math.sqrt(2))*np.ones((n_pts))
    edgy           = 1.0*np.ones((n_pts))
    ctrl_pts       = generate_cylinder_pts(n_pts)
    ctrl_pts[:,:] *= shape_size

if (shape_type == 'square'):
    radius         = np.zeros((n_pts))
    edgy           = np.ones((n_pts))
    ctrl_pts       = generate_square_pts(n_pts)
    ctrl_pts[:,:] *= shape_size

if (shape_type == 'random'):
    radius         = np.random.uniform(low=0.0, high=1.0, size=n_pts)
    edgy           = np.random.uniform(low=0.0, high=1.0, size=n_pts)
    ctrl_pts       = np.random.rand(n_pts,2)
    ctrl_pts[:,:] *= shape_size

shape = Shape(shape_name,
              ctrl_pts,
              n_pts,
              n_sampling_pts,
              radius,
              edgy,
              lattice.output_dir)

shape.generate()
shape.generate_image()
shape.write_csv()

# Generate lattice from shape closed curve
lattice.add_obstacle(shape.curve_pts)
lattice.generate_image()

# Initialize fields
lattice.set_full_poiseuille(u_lbm, rho_lbm)
lattice.u[:,lattice.lattice] = 0.0

# Set initial distributions
lattice.equilibrium()
lattice.g = lattice.g_eq.copy()

# Solve
bar = progress.bar.Bar('Solving...', max=it_max)
for it in range(it_max+1):

    # Compute macroscopic fields
    lattice.macro()

    # Output field
    lattice.output_fields(it,
                          output_freq,
                          u_norm   = True,
                          u_stream = False)

    # Compute equilibrium state
    lattice.equilibrium()

    # Collisions
    lattice.trt_collisions()

    # Streaming
    lattice.stream()

    # Boundary conditions
    lattice.bounce_back_obstacle()
    lattice.zou_he_bottom_wall_velocity()
    lattice.zou_he_left_wall_velocity()
    lattice.zou_he_right_wall_pressure()
    lattice.zou_he_top_wall_velocity()
    lattice.zou_he_bottom_left_corner()
    lattice.zou_he_top_left_corner()
    lattice.zou_he_top_right_corner()
    lattice.zou_he_bottom_right_corner()

    # Compute drag/lift
    lattice.drag_lift(it, rho_lbm, u_avg, D_lbm)

    # Increment bar
    bar.next()

# End bar
bar.finish()
