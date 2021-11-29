# Generic imports
import os
import pytest

# Custom imports
from lbm.src.app.turek    import *
from lbm.src.core.lattice import *

###############################################
### Test lattice obstacle generation
def test_lattice():

    # Initial space
    print("")

    # Generate turek app and lattice
    app = turek()
    app.L_lbm = 100
    app.compute_lbm_parameters()
    ltc = lattice(app)

    # Add obstacle
    app.add_obstacles(ltc, app.obstacles)
    n_bnd = app.obstacles[0].boundary.shape[0]

    # Check nb of boundary pts
    assert(n_bnd == 234)

    # Increase accuracy and redo
    app = turek()
    app.L_lbm = 200
    app.compute_lbm_parameters()
    ltc = lattice(app)
    app.add_obstacles(ltc, app.obstacles)
    n_bnd = app.obstacles[0].boundary.shape[0]
    assert(n_bnd == 468)
