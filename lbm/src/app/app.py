# Custom imports
from lbm.src.core.factory import *
from lbm.src.app.cavity   import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("cavity", cavity)
