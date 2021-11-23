# Custom imports
from lbm.src.core.factory import *
from lbm.src.app.cavity   import *
from lbm.src.app.turek    import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("cavity", cavity)
app_factory.register("turek",  turek)
