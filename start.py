# Generic imports
import os
import sys
import time
import numpy as np

# Custom imports
from lbm.src.app.app import *

########################
# Run lbm simulation
########################
if __name__ == '__main__':

    # Check command-line input
    if (len(sys.argv) == 2):
        app_name = sys.argv[1]
    else:
        print('Command line error, please use as follows:')
        print('python3 start.py app_name')

    # Instanciate app
    app = app_factory.create(app_name)

    # Instanciate lattice
    lattice = lattice(app)

    # Initialize fields and distributions
    lattice.initialize()
    lattice.equilibrium()
    lattice.g = lattice.g_eq.copy()

    # Count time
    start_time = time.time()

    # Solve
    it      = 0
    compute = True
    print('### Solving')
    while (compute):

        # Printings
        lattice.printings(it)

        # Set inlets
        lattice.set_inlets()

        # Compute macroscopic fields
        lattice.macro()

        # Output field
        lattice.outputs(it)

        # Compute equilibrium state
        lattice.equilibrium()

        # Streaming
        lattice.collision_stream()

        # Boundary conditions
        lattice.set_bc()

        # Compute observables (drag, lift, etc)
        lattice.observables()

        # Check stopping criterion
        compute = lattice.check_stop(it)

        # Increment iteration
        it += 1

# Count time
end_time = time.time()
print("# Loop time = {:f}".format(end_time - start_time))

# Perform final operations and outputs
lattice.finalize()
