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
    app.initialize(lattice)

    # Timer and loop data
    start_time = time.time()
    it         = 0
    compute    = True

    # Solve
    print('### Solving')
    while (compute):

        # Printings
        app.printings(it)

        # Set inlets
        app.set_inlets(lattice, it)

        # Compute macroscopic fields
        lattice.macro()

        # Output field
        app.outputs(lattice, it)

        # Compute equilibrium state
        lattice.equilibrium()

        # Streaming
        lattice.collision_stream()

        # Boundary conditions
        app.set_bc(lattice)

        # Compute observables (drag, lift, etc)
        app.observables(lattice, it)

        # Check stopping criterion
        compute = app.check_stop(it)

        # Increment iteration
        it += 1

# Count time
end_time = time.time()
print("# Loop time = {:f}".format(end_time - start_time))

# Perform final operations and outputs
app.finalize(lattice)
