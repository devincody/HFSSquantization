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
    nValues = 4
    values = np.linspace(0.0,.0003,nValues)
    scale_z = 0.000
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
        sim_wg.build_L_mat(verbose = False)
        sim_wg.build_C_mat(verbose = False)
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
        diff = np.abs(100*(freq-modes)/modes)
        print "Difference:", diff
        frequency_vector.append(diff)
        
    hfss.release() 
    print "optimize scalex"
    print values
    print frequency_vector
    np.save("../data/frequencyvector", frequency_vector)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(values, frequency_vector)
    plt.xlabel("scalex")
    plt.ylabel("error")  
    plt.show()
    
def optimize_scalez():
    print "optimize scalez"
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    nValues = 3
    values = np.linspace(0.000,.0001,nValues)
    scale_x = 0.0001
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
        diff = 100*(freq-modes)/modes
        print "Difference:", diff
        frequency_vector.append(diff)
        
    hfss.release() 
    print "optimize scalez"
    print values
    print frequency_vector
    np.save("../data/frequencyvector", frequency_vector)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(121)
    ax1.plot(values, frequency_vector)
    ax2 = fig1.add_subplot(122)
    ax2.plot(values, np.abs(frequency_vector))
    plt.xlabel("scalez")
    plt.ylabel("error")  
    plt.show()
    
    
def One_Trial(N = 100, collect_new_data = True, qubit = False, qu_theta = 0, verbose = False):
    print "One Trial"
    Load_sim_number = 1
    wg = waveguide() 
    frequency_vector = []
    scale_x = 0.0001
    wg.set_scalex(scale_x)
    z = 0.0001
    wg.set_scalez(z)
    wg = waveguide(angle_s = 0, angle_e = 360,angle_n = N) 
    if collect_new_data:   
        wg.compute_LCVI(cap_surf = "CrossSecIntSurf", ind_surf = "CrossSecIntSurf1")  
        wg.save(Load_sim_number)
    else:
        wg.load(Load_sim_number)

    #wg.plot()
    wg.update_HFSS_Solutions()
    modes = wg.eigenmodes[0][0:2]
    sim_wg = simulated_wg(Load_sim_number)

    #GET LC VALUES and BUILD MATRIX
    sim_wg.build_L_mat(N, qubit, qu_theta,verbose)
    sim_wg.build_C_mat_parallel(N, qubit, qu_theta,verbose)
    freq_tmp = sim_wg.get_frequencies()
    hfss.release()
    #Get best values
    i = 0
    while (freq_tmp[i] < 5E9):
        i++
        
    freq = freq_tmp[i:i+2]/10**9
    print "%d:%d"%(i,i+2)
        
    #Display Results
    
    print "scalex", scale_x
    print "scalez (opt)", z
    if verbose:
        print freq_tmp
        print "Eigenmode Frequencies:", freq


    print "HFSS Frequencies:", modes
    diff = 100*(freq-modes)/modes
    print "Percent Difference:", diff
    frequency_vector.append(diff)
    np.save("../data/frequencyvector", freq_tmp)
    return freq,modes
   
def multi_perturb():
    N =100
    #One_Trial(100)
    freq= []
    for z_angle in range(N):        
        f = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle)#(360-z_angle)/360*N)
        print f        
        freq.append(f)
    
    plt.plot(range(N),freq)
    plt.show()    
    
def collect_and_perturb():
    N =100
    #One_Trial(100)
    R = 10
    z_angles = np.linspace(0,N-R/N,R)
    freq_eigen= []
    freq_hfss = []
    for z_angle in z_angles:       
        design = hfss.get_active_design()
        design.set_variable('zmonAngle',(u'%.4f' % (z_angle)))
        f,m = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle,verbose = False)#(360-z_angle)/360*N)
        print f,m       
        freq_eigen.append(f)
        freq_hfss.append(m)
    
    print "Z_angles", z_angles
    print "freq eigen", freq_eigen
    print "freq hfss", freq_hfss
    
    plt.plot(z_angles,freq_eigen)
    plt.show()    
    plt.plot(z_angles,freq_hfss)
    plt.show()    
    plt.plt(z_angles,freq_eigen,z_angles,freq_hfss)
    plt.show()
    
def single_perturb():
    design = hfss.get_active_design()
    z_angle = float(design.get_variable_value('zmonAngle'))
    N =100
    freq= []
    One_Trial(N, collect_new_data = False, qubit = True, qu_theta = (360-z_angle)/360*N,verbose = True)

    print "Zmon_angle", z_angle

if __name__ == "__main__":
    #One_Trial(N, collect_new_data = True, qubit = False, qu_theta = 0*z_angle,verbose = False)
    collect_and_perturb()