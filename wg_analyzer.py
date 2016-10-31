# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:23:19 2016
@author: Devin

Code that Utilizes functions in wg_eigenmode_simulator.py and LCVI_array_builder.py
to analyze some properties of the waveguide
"""

import sys
if ~(r'F:\Documents\Yale\Devoret_Research\HFSSpython\pyHFSS' in sys.path):
    sys.path.append(r'F:\Documents\Yale\Devoret_Research\HFSSpython\pyHFSS')
import hfss
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from wg_eigenmode_simulator import simulated_wg
from LCVI_array_builder import waveguide

def main():
    Load_sim_number = 1
    wg = waveguide(angle_n = 20)   
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
    print "optimize scalex"
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    nValues = 11
    values = np.linspace(0.0,.001,nValues)
    scale_z = 0.0005
    wg.set_scalez(scale_z)
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
        sim_wg.build_L_mat(verbose = True)
        sim_wg.build_C_mat(verbose = True)
        freq_tmp = sim_wg.get_frequencies()
        if (freq_tmp[0] > 3*(10**9)):
            freq = freq_tmp[0:2]/10**9
            print freq_tmp[0]
            print "0:2"
        else:
            freq = freq_tmp[1:3]/10**9
            print"1:3"
        print "scalex (opt)", x
        print "scalez", scale_z
        print "Eigenmode Frequencies:", freq
        print "HFSS Frequencies:", modes
        diff = (freq-modes)/freq
        print "Difference:", diff
        frequency_vector.append(diff)
        
    hfss.release() 
    print "optimize scalex"
    print values
    print frequency_vector
    np.save("../data/frequencyvector", frequency_vector)
    plt.plot(values, frequency_vector)
    plt.show()
    
def optimize_scalez():
    print "optimize scalez"
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    nValues = 11
    values = np.linspace(0.0000,.0005,nValues)
    scale_x = 0.0005
    wg.set_scalex(scale_x)
    for i in range(nValues):
        wg = waveguide() 
        z = values[i]
        wg.set_scalez(z)
        wg.compute_LCVI()  
        #wg.load(Load_sim_number)
        wg.plot()
        wg.save(Load_sim_number)
        #
        #wg.plot()
        modes = wg.eigenmodes[0][0:2]
        
        sim_wg = simulated_wg(Load_sim_number)
        #sim_wg.test_interpolate()
        sim_wg.build_L_mat(verbose = True)
        sim_wg.build_C_mat(verbose = True)
        freq_tmp = sim_wg.get_frequencies()
        if (freq_tmp[0] > 3*(10**9)):
            freq = freq_tmp[0:2]/10**9
            print freq_tmp[0]
            print "0:2"
        else:
            freq = freq_tmp[1:3]/10**9
            print"1:3"
        #print freq_tmp
        print "scalex", scale_x
        print "scalez (opt)", z
        print "Eigenmode Frequencies:", freq
        print "HFSS Frequencies:", modes
        diff = (freq-modes)/freq
        print "Difference:", diff
        frequency_vector.append(diff)
        
    hfss.release() 
    print "optimize scalez"
    print values
    print frequency_vector
    np.save("../data/frequencyvector", frequency_vector)
    plt.plot(values, frequency_vector)
    plt.show()
    
def One_Trial():
    print "One Trial"
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    scale_x = 0.0005
    wg.set_scalex(scale_x)
    z = 0.00004
    wg.set_scalez(z)
    wg = waveguide(angle_n = 100) 
    wg.compute_LCVI(cap_surf = "CrossSecIntSurf", ind_surf = "CrossSecIntSurf1")  
    #wg.load(Load_sim_number)
    #wg.plot()
    wg.save(Load_sim_number)
    #
    #wg.plot()
    modes = wg.eigenmodes[0][0:2]
    
    sim_wg = simulated_wg(Load_sim_number)
    #sim_wg.test_interpolate()
    
    #GET LC VALUES and BUILD MATRIX
    sim_wg.build_L_mat(verbose = True)
    sim_wg.build_C_mat(verbose = True)
    freq_tmp = sim_wg.get_frequencies()
    hfss.release() 
    #Get best values
    if (freq_tmp[0] > 3*(10**9)):
        freq = freq_tmp[0:2]/10**9
        print freq_tmp[0]
        print "0:2"
    else:
        freq = freq_tmp[1:3]/10**9
        print"1:3"
        
    #Display Results
    print freq_tmp # print all HFSS Frequencies
    print "scalex", scale_x
    print "scalez (opt)", z
    print "Eigenmode Frequencies:", freq
    print "HFSS Frequencies:", modes
    diff = (freq-modes)/freq
    print "Difference:", diff
    frequency_vector.append(diff)
        
    np.save("../data/frequencyvector", freq_tmp)
   
    
    
if __name__ == "__main__":
    #One_Trial()
    optimize_scalex()
    #optimize_scalez()
    #main()    