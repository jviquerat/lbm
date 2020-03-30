# Generic imports
import os
import math
import numpy as np

###############################################
### Class Buff
### A buffer class
class Buff:

    ### Create object
    def __init__(self,
                 name,
                 dt,
                 output_dir):

        # Fill structure
        self.name         = name
        self.buff         = np.zeros([2])
        self.avg1_buff    = np.zeros([2])
        self.avg2_buff    = np.zeros([2])
        self.avg3_buff    = np.zeros([2])
        self.it           = 0
        self.dt           = dt
        self.output_dir   = output_dir
        self.freq         = 1.0
        self.obs          = 0.0
        self.obs_cv_cnt   = 0
        self.obs_cv       = False

    ### Add a value to the buffer
    def add(self, value):

        self.buff = np.append(self.buff, value)
        self.it  += 1

    ### Full average buffer
    def f_avg(self):

        return np.sum(self.buff)/len(self.buff)

    ### Partial average buffer
    def p_avg(self, buff, i ,j):

        return np.sum(buff[i:j])/(float(j-i+1))

    ### Compute average of moving average
    def mv_avg(self):

        f_avg             = self.f_avg()
        it_s              = math.floor(2*self.it/3)
        it_e              = self.it
        self.obs          = self.p_avg(self.buff, it_s, it_e)
        self.avg1_buff    = np.append(self.avg1_buff, self.obs)
        self.obs          = self.p_avg(self.avg1_buff, it_s, it_e)
        self.avg2_buff    = np.append(self.avg2_buff, self.obs)
        self.obs          = self.p_avg(self.avg2_buff, it_s, it_e)
        self.avg3_buff    = np.append(self.avg3_buff, self.obs)
        self.obs          = self.p_avg(self.avg3_buff, it_s, it_e)

        growth = 0.0
        if (self.it > 4):
            growth = (self.avg3_buff[-1] -
                      self.avg3_buff[-4])/(3.0*self.dt)

            if (abs(growth) < 1.0e-2):
                self.obs_cv_cnt += 1
            else:
                self.obs_cv_cnt  = 0

            if (self.obs_cv_cnt > 1000):
                self.obs_cv = True

        return self.obs, growth
