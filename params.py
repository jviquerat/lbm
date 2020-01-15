# Generic imports
from datetime import datetime

###############################################
### Set CFD/DRL parameters
###############################################

# Set CFD parameters
xmin           =-10.0             # for domain size
xmax           = 20.0             # for domain size
ymin           =-5.0              # for domain size
ymax           = 5.0              # for domain size
u_in           = 1.0              # input velocity
rho            = 1.0              # fluid density
timestep       = 0.5              # cfd timestep
cfl            = 0.25             # CFL
reynolds       = 100.0            # Re value for L=1
domain_type    = 'channel'        # 'open' or 'channel'
name           = 'shape'          # descriptive name for outputs

# Set general parameters
n_cpu          = 8               # nb of cpus
reset_dir      = 'reset/4/'       # reset dir for shape
pts_to_move    = [0]        # define moving pts
n_params       = len(pts_to_move) # nb of moving pts
min_rad        = 0.1              # min radius
max_rad        = 2.0              # max radius
min_edg        = 0.0              # min edgy
max_edg        = 1.0              # max edgy
shape_h        = 0.5              # mesh accuracy around shape
domain_h       = 0.5              # mesh accuracy in domain
n_best         = n_cpu            # nb of best shapes

# Set PPO parameters
max_generations = 50              # max nb of episodes
policy_arch     = [2,2]           # architecture of agent
learning_freq   = 1               # update frequency for agent
learning_rate   = 5.0e-3          # learning rate for policy net
clip_range      = 0.5             # clipping value for PPO
gamma           = 0.9             # discount factor
rwd_shaping     = 5.0             # Shaping value for reward
rwd_penal       =-5.0             # Max penalty for reward
n_mini_batches  = 1               # minibatch size
n_epochs        = 32 #32 ou 64 max    # nb of epochs for the drl net

# Set output parameters
time          = datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
output_path   = 'outputs/'+name+'_'+time+'/'
mesh_path     = 'mesh/'
drag_path     = 'drag/'
sol_path      = 'sol/'
reward_path   = './'
action_path   = './'
state_path    = './'
csv_path      = 'csv/'
shape_path    = 'shape/'
cfd_verbose   = False      # CFD verbosity
make_video    = False      # output vtu at each time-step

# Set sorting parameters
sorted_dir = 'sorted/'+name+'_'+time+'/'
sol_dir    = sorted_dir+'sol/'
best_dir   = sorted_dir+'best/'
csv_dir    = sorted_dir+'csv/'
shape_dir  = sorted_dir+'shape/'
ring_size  = 100
