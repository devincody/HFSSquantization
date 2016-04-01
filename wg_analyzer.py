# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:23:19 2016
@author: Devin

Code that Utilizes functions in wg_eigenmode_simulator.py adn LCVI_array_builder.py
to analyze some properties of the waveguide
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

def main():
    Load_sim_number = 1
    wg = waveguide()   
    wg.compute_LCVI()  
    wg.save(Load_sim_number)
    wg.load(Load_sim_number)
    wg.plot()
    modes = wg.eigenmodes[0][0:2]
    
    sim_wg = simulated_wg(Load_sim_number)
    sim_wg.test_interpolate()
    sim_wg.build_L_mat(verbose=True)
    sim_wg.build_C_mat(verbose=True)
    freq = sim_wg.get_frequencies()[1:3]/10**9
    print "Simulated Frequencies:", freq
    print "HFSS Frequencies:", modes
    diff = (freq-modes)/freq
    print "Difference:", diff
        
    hfss.release() 

def optimize_scalex():
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    nValues = 11
    values = np.linspace(0.0,.001,nValues)
    for i in range(nValues):
        wg = waveguide() 
        x = values[i]
        wg.set_scalex(x)
        wg.compute_LCVI()      
        wg.save(Load_sim_number)
        #wg.load(Load_sim_number)
        #wg.plot()
        modes = wg.eigenmodes[0][0:2]
        
        sim_wg = simulated_wg(Load_sim_number)
        #sim_wg.test_interpolate()
        sim_wg.build_L_mat()
        sim_wg.build_C_mat()
        freq = sim_wg.get_frequencies()[1:3]/10**9
        print "Simulated Frequencies:", freq
        print "HFSS Frequencies:", modes
        diff = (freq-modes)/freq
        print "Difference:", diff
        frequency_vector.append(freq)
        
    hfss.release() 
    print frequency_vector
    np.save("../data/frequencyvector", frequency_vector)
    plt.plot(values, frequency_vector)
    plt.show()
    

if __name__ == "__main__":
    main()    