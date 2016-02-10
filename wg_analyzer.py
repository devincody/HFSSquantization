# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:23:19 2016
@author: Devin
"""

import sys
if ~(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS' in sys.path):
    sys.path.append(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS')
import hfss, numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from wg_simulator import simulated_wg
from parametricTest import waveguide

wg = waveguide() 
#wg.compute_LCVI()      
#print wg.inductance
#wg.save()
wg.load()   
wg.plot()
modes = wg.eigenmodes[0]
hfss.release()
sim_wg = simulated_wg()
#sim_wg.test_interpolate()
sim_wg.build_L_mat()
sim_wg.build_C_mat()
freq = sim_wg.get_frequencies()[1:3]/10**9
print "Simulated Frequencies:", freq
print "HFSS Frequencies:", modes
diff = (freq-modes)/freq
print diff