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

