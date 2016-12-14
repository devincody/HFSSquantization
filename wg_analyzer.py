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
    n_modes = 3
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
    modes = wg.eigenmodes[0][0:n_modes]
    sim_wg = simulated_wg(Load_sim_number)

    #GET LC VALUES and BUILD MATRIX
    sim_wg.build_L_mat(N, qubit, qu_theta,verbose)
    sim_wg.build_C_mat_parallel(N, qubit, qu_theta,verbose)
    freq_tmp = sim_wg.get_frequencies()
    hfss.release()
    #Get best values
#    i = 0
#    while (freq_tmp[i] < 7E9):
#        i =i+1
    if freq_tmp[0]> 7.1E9:
        i = 0
        freq2 = freq_tmp[i:i+2]/10**9
        freq2.insert(0,'nan')
    
    i = 0
    freq = freq_tmp[i:i+n_modes]/10**9
    print "%d:%d"%(i,i+n_modes)
        
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
    return freq,modes,diff, freq2
   
def multi_perturb():
    N =100
    #One_Trial(100)
    freq_eigen= []
    diff = []
    for z_angle in range(N):        
        f,m,d = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle)#(360-z_angle)/360*N)
        print f,m,d   
        freq_eigen.append(f)

    z_angles = np.linspace(0,99,100)
    freq_hfss = [[7.17949645918, 7.52311254486], [7.17729060558, 7.52642053675], [7.18020417597, 7.52492349072], [7.18751505185, 7.51907559426], [7.19484172028, 7.51965879766], [7.2018582451, 7.51362016454], [7.21203045676, 7.51426117744], [7.22111698839, 7.50675208528], [7.232559786, 7.49761448872], [7.24373263931, 7.49354303057], [7.25169405444, 7.48509848577], [7.26153969295, 7.48479273155], [7.2685503918, 7.4779029865], [7.27666774887, 7.47384766439], [7.28592255107, 7.46257280371], [7.2884733008, 7.46779290153], [7.29147749386, 7.45915859387], [7.29553744083, 7.45928289], [7.29372978243, 7.4568602345], [7.29956383406, 7.46160199631], [7.30250987478, 7.46128985878], [7.29984683185, 7.46222505765], [7.30177099813, 7.46834402144], [7.30009053036, 7.47076497162], [7.29772867153, 7.4673487338], [7.29706466923, 7.46970499098], [7.2931305053, 7.47621135583], [7.29456884529, 7.48096259899], [7.29188152292, 7.48274498598], [7.29241676633, 7.48854175253], [7.28784841727, 7.49298119333], [7.28756283979, 7.49547526009], [7.28164647066, 7.49587901912], [7.28025619732, 7.49883869637], [7.28050585401, 7.50296030417], [7.27760409152, 7.50752801602], [7.27444670945, 7.51282189643], [7.27059386037, 7.51390330473], [7.27298649252, 7.51419114235], [7.26602532856, 7.51264916413], [7.269937045, 7.51919104958], [7.26548863661, 7.52045287139], [7.26342014554, 7.52113092101], [7.26666251186, 7.52299068981], [7.26625892072, 7.52129119286], [7.26469938545, 7.52333746679], [7.26238470764, 7.51730940131], [7.26197109232, 7.51705621858], [7.26391515593, 7.52203085334], [7.26149144728, 7.51819853355], [7.25970003983, 7.5206041302], [7.26056233098, 7.52111598943], [7.26218407255, 7.52119912248], [7.26183550999, 7.52080346143], [7.26070655725, 7.5207470539], [7.26385010518, 7.5271399089], [7.26530066781, 7.52033012405], [7.26440787125, 7.51566960356], [7.26755512969, 7.51979501886], [7.26716562431, 7.51685406984], [7.27241655815, 7.51611041641], [7.27097404378, 7.51046579307], [7.26993409032, 7.51113765311], [7.27519148604, 7.51231760159], [7.27320133048, 7.50736769254], [7.27845028092, 7.49865766676], [7.27494047469, 7.49602312553], [7.2834117512, 7.49695314167], [7.28218697902, 7.49346843524], [7.28689137486, 7.49690334089], [7.28614937781, 7.48627108648], [7.29218377844, 7.4891743447], [7.2948763783, 7.48240965356], [7.29513594687, 7.48333279004], [7.29696093364, 7.48109753616], [7.2973703338, 7.47181520478], [7.29903960354, 7.47533198758], [7.29902339551, 7.4669961955], [7.30014214629, 7.46357673546], [7.30273549286, 7.46120629625], [7.29979030443, 7.46186697785], [7.30266969891, 7.46163653316], [7.29670037061, 7.46225370245], [7.29679233747, 7.46221727325], [7.28897571951, 7.46070673436], [7.28886570178, 7.46822920942], [7.28411671113, 7.46896567754], [7.27456017031, 7.4713237753], [7.2675392359, 7.48177948259], [7.25926750472, 7.48394492336], [7.24408352307, 7.48496967572], [7.23721494457, 7.49086989952], [7.23273524807, 7.50057244084], [7.21753166348, 7.50757900796], [7.21196444948, 7.50993038214], [7.19846553379, 7.51245722219], [7.194024376, 7.5163601517], [7.18289267046, 7.5256934201], [7.18400433112, 7.52457543177], [7.17566000817, 7.53041260336]]    
    
    for i in range(N):
        diff.append([100*(freq_eigen[i][0] - freq_hfss[i][0])/freq_hfss[i][0], 100*(freq_eigen[i][1] - freq_hfss[i][1])/freq_hfss[i][1]])
    
    plt.plot((z_angles*360/N),freq_eigen)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Eigenmode Frequencies [GHz]")
    plt.title("Perturbatively Calculated Eigenmodes")
    plt.show()    
    plt.plot((z_angles*360/N),freq_hfss)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("HFSS Frequencies [GHz]")
    plt.title("HFSS Calculated Eigenmodes")
    plt.show()    
    plt.plot((z_angles*360/N),freq_eigen,'r',(z_angles*360/N),freq_hfss,'b')
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Frequencies [GHz]")
    plt.title("Perturbatively vs. HFSS Eigenmodes")
    plt.show()
    plt.plot((z_angles*360/N),diff)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Difference in Frequencies [Percent]")
    plt.title("Perturbatively vs. HFSS Eigenmodes ERROR")
    plt.show()  

    print sum(np.abs(diff))/len(diff)
    
def collect_and_perturb():
    N =100
    #One_Trial(100)
    R = 11
    z_angles = np.linspace(0,50,R)
    freq_eigen= []
    freq_hfss = []
    diff = []
    for z_angle in z_angles:       
        design = hfss.get_active_design()
        design.set_variable('zmonAngle',(u'%.4f' % (z_angle*360/N)))
        f,m,d = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle,verbose = False)#(360-z_angle)/360*N)
        print f,m,d   
        freq_eigen.append(f)
        freq_hfss.append(m)
        diff.append(d)
    
    print "Z_angles", z_angles
    print "freq eigen", freq_eigen
    print "freq hfss", freq_hfss
    print "diff", diff
    
    plt.plot((z_angles*360/N),freq_eigen)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Eigenmode Frequencies [GHz]")
    plt.title("Perturbatively Calculated Eigenmodes")
    plt.show()    
    plt.plot((z_angles*360/N),freq_hfss)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("HFSS Frequencies [GHz]")
    plt.title("HFSS Calculated Eigenmodes")
    plt.show()    
    plt.plot((z_angles*360/N),freq_eigen,'r',(z_angles*360/N),freq_hfss,'b')
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Frequencies [GHz]")
    plt.title("Perturbatively vs. HFSS Eigenmodes")
    plt.show()
    plt.plot((z_angles*360/N),diff)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Difference in Frequencies [GHz]")
    plt.title("Perturbatively vs. HFSS Eigenmodes ERROR")
    plt.show()
    
def single_perturb():
    design = hfss.get_active_design()
    z_angle = float(design.get_variable_value('zmonAngle'))
    N =100
    freq= []
    One_Trial(N, collect_new_data = False, qubit = True, qu_theta = (360-z_angle)/360*N,verbose = True)

    print "Zmon_angle", z_angle

if __name__ == "__main__":
    #One_Trial(N, collect_new_data = False, qubit = True, qu_theta = 0*z_angle,verbose = True)
    #collect_and_perturb()
    multi_perturb()