# Generic imports
import os
import pytest

# Custom imports
from lbm.src.app.poiseuille import *
from lbm.src.core.lattice   import *
from lbm.src.core.run       import *

###############################################
### Test poiseuille
def test_poiseuille():

    # Initial space
    print("")

    # Run simple poiseuille flow
    app = poiseuille()
    ltc = lattice(app)
    run(ltc, app)

    # Check error
    l1_error = app.compute_error(ltc)
    print(l1_error)
    assert(l1_error < 1.0e-3)
