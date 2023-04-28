# Generic imports
import os
import pytest

# Custom imports
from lbm.src.app.cavity   import *
from lbm.src.core.lattice import *
from lbm.src.core.run     import *

###############################################
### Test cavity
def test_cavity():

    # Initial space
    print("")

    # Run simple driven cavity
    app = cavity()
    app.output_freq = 10000
    ltc = lattice(app)
    run(ltc, app)

    # Check error
    vx, uy = app.line_fields(ltc)
    assert(abs(vx[10] - 0.12493089684236539) < 1.0e-6)
    assert(abs(vx[50] - 0.05295104908939561) < 1.0e-6)
    assert(abs(uy[10] + 0.05968571489630510) < 1.0e-6)
    assert(abs(uy[50] + 0.19792323493599165) < 1.0e-6)
