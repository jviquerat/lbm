# Generic imports
import os
import math
import numpy as np

###############################################
### Class Buffer
### A buffer class
class Buffer:

    ### Create object
    def __init__(self,
                 size):

        # Check inputs
        if (size < 1):
            raise ValueError('Incorrect buffer size')

        if (size is not int):
            size = math.floor(size)

        # Fill structure
        self.size  = size
        self.buff  = np.zeros([self.size])
        self.index = 0

    ### Add a value to the buffer
    def add(self,
            value):

        # Add value and increment counter
        self.buff[self.index] = value
        self.index           += 1

        # Circular test
        if (self.index == self.size): self.index = 0

    ### Average buffer
    def avg(self):

        # Average and return
        return np.sum(self.buff)/self.size
