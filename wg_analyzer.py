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
    '''
    Code builds LC Matricies and calculates eigenfrequencies. Optionally includes qubit and can also 
    Resimulate WGMR for HFSS frequencies, and recalculate LCVI data as needed.
    
    This function is called multiple times by multi_perturb to obtain a complete picture of the 
    quantum parameters as a function of qubit angle    
    '''
    print "One Trial"
    Load_sim_number = 1
    n_modes = 5
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
        wg.update_HFSS_Solutions()
    else:
        wg.load(Load_sim_number)

    #wg.plot()
    
    modes = [0,0,0,0,0]#wg.eigenmodes[0][0:n_modes]
    sim_wg = simulated_wg(Load_sim_number)

    #GET LC VALUES and BUILD MATRIX
    sim_wg.build_L_mat(N, qubit, qu_theta, verbose)
    sim_wg.build_C_mat_parallel(N, qubit, qu_theta, verbose)
    freq_tmp,chi = sim_wg.get_frequencies(qu_theta)
    hfss.release()
    
    #Get best values
#    i = 0
#    while (freq_tmp[i] < 7E9):
#        i =i+1
    i = 0
    freq = freq_tmp[i:i+n_modes]/10**9
    if freq_tmp[0]< 4E9:
        i=1
        freq2 = list(freq_tmp[i:i+3]/10**9)
    else:
        freq2 = freq
    #print "%d:%d"%(i,i+n_modes)
    
    diff = 100*(freq-modes)/modes    
    frequency_vector.append(diff)
    #Display Results        
    if verbose:
        print freq_tmp
        print "Eigenmode Frequencies:", freq
        print "HFSS Frequencies:", modes
        print "Percent Difference:", diff
        
    np.save("../data/frequencyvector", freq_tmp)
    return freq,modes,diff, freq2,chi
   
def multi_perturb():
    N = 100 #number of angles around WGMR
    #One_Trial(100)
    freq_eigen= [] #finding the eigen frequencies via LC Matrix
    diff = []
    alpha = []
    chi_2 = []
    for z_angle in range(N):        
        f,m,d,f2,chi = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle)#(360-z_angle)/360*N)
        print f,m,d   
        freq_eigen.append(f2[1:3])
        alpha.append([chi[0][0]/1E6,chi[1][1]/1E6,chi[2][2]/1E6])
        chi_2.append([chi[0][1]/1E6,chi[0][2]/1E6])

    z_angles = np.linspace(0,99,N)
    freq_hfss = [[7.17949645918, 7.52311254486], [7.17729060558, 7.52642053675], [7.18020417597, 7.52492349072], [7.18751505185, 7.51907559426], [7.19484172028, 7.51965879766], [7.2018582451, 7.51362016454], [7.21203045676, 7.51426117744], [7.22111698839, 7.50675208528], [7.232559786, 7.49761448872], [7.24373263931, 7.49354303057], [7.25169405444, 7.48509848577], [7.26153969295, 7.48479273155], [7.2685503918, 7.4779029865], [7.27666774887, 7.47384766439], [7.28592255107, 7.46257280371], [7.2884733008, 7.46779290153], [7.29147749386, 7.45915859387], [7.29553744083, 7.45928289], [7.29372978243, 7.4568602345], [7.29956383406, 7.46160199631], [7.30250987478, 7.46128985878], [7.29984683185, 7.46222505765], [7.30177099813, 7.46834402144], [7.30009053036, 7.47076497162], [7.29772867153, 7.4673487338], [7.29706466923, 7.46970499098], [7.2931305053, 7.47621135583], [7.29456884529, 7.48096259899], [7.29188152292, 7.48274498598], [7.29241676633, 7.48854175253], [7.28784841727, 7.49298119333], [7.28756283979, 7.49547526009], [7.28164647066, 7.49587901912], [7.28025619732, 7.49883869637], [7.28050585401, 7.50296030417], [7.27760409152, 7.50752801602], [7.27444670945, 7.51282189643], [7.27059386037, 7.51390330473], [7.27298649252, 7.51419114235], [7.26602532856, 7.51264916413], [7.269937045, 7.51919104958], [7.26548863661, 7.52045287139], [7.26342014554, 7.52113092101], [7.26666251186, 7.52299068981], [7.26625892072, 7.52129119286], [7.26469938545, 7.52333746679], [7.26238470764, 7.51730940131], [7.26197109232, 7.51705621858], [7.26391515593, 7.52203085334], [7.26149144728, 7.51819853355], [7.25970003983, 7.5206041302], [7.26056233098, 7.52111598943], [7.26218407255, 7.52119912248], [7.26183550999, 7.52080346143], [7.26070655725, 7.5207470539], [7.26385010518, 7.5271399089], [7.26530066781, 7.52033012405], [7.26440787125, 7.51566960356], [7.26755512969, 7.51979501886], [7.26716562431, 7.51685406984], [7.27241655815, 7.51611041641], [7.27097404378, 7.51046579307], [7.26993409032, 7.51113765311], [7.27519148604, 7.51231760159], [7.27320133048, 7.50736769254], [7.27845028092, 7.49865766676], [7.27494047469, 7.49602312553], [7.2834117512, 7.49695314167], [7.28218697902, 7.49346843524], [7.28689137486, 7.49690334089], [7.28614937781, 7.48627108648], [7.29218377844, 7.4891743447], [7.2948763783, 7.48240965356], [7.29513594687, 7.48333279004], [7.29696093364, 7.48109753616], [7.2973703338, 7.47181520478], [7.29903960354, 7.47533198758], [7.29902339551, 7.4669961955], [7.30014214629, 7.46357673546], [7.30273549286, 7.46120629625], [7.29979030443, 7.46186697785], [7.30266969891, 7.46163653316], [7.29670037061, 7.46225370245], [7.29679233747, 7.46221727325], [7.28897571951, 7.46070673436], [7.28886570178, 7.46822920942], [7.28411671113, 7.46896567754], [7.27456017031, 7.4713237753], [7.2675392359, 7.48177948259], [7.25926750472, 7.48394492336], [7.24408352307, 7.48496967572], [7.23721494457, 7.49086989952], [7.23273524807, 7.50057244084], [7.21753166348, 7.50757900796], [7.21196444948, 7.50993038214], [7.19846553379, 7.51245722219], [7.194024376, 7.5163601517], [7.18289267046, 7.5256934201], [7.18400433112, 7.52457543177], [7.17566000817, 7.53041260336]]    
    angles_hfss = np.linspace(0,99,100) *360/100   
    
    for i in range(N):
        diff.append([100*(freq_eigen[i][0] - freq_hfss[i][0])/freq_hfss[i][0], 100*(freq_eigen[i][1] - freq_hfss[i][1])/freq_hfss[i][1]])
          
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('linear')
    ax.plot((z_angles*360/N),freq_eigen)

    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("Eigenmode Frequencies [GHz]")
    ax.set_xlim([0,360])
    ax.set_title("Perturbatively Calculated Eigenmodes")
    fig.show()  

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('linear')
    ax.plot(angles_hfss,freq_hfss)
    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("HFSS Frequencies [GHz]")
    ax.set_xlim([0,360])
    ax.set_title("HFSS Calculated Eigenmodes")
    fig.show()  
    


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('linear')
    ax.plot((z_angles*360/N),freq_eigen,'r',angles_hfss,freq_hfss,'b')
    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("Frequencies [GHz]")
    ax.set_xlim([0,360])
    ax.set_title("Perturbatively vs. HFSS Eigenmodes")
    fig.show()  
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('linear')
    ax.plot((z_angles*360/N),diff)
    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("Difference in Frequencies [%]")
    ax.set_xlim([0,360])
    ax.set_title("Perturbatively vs. HFSS Eigenmodes ERROR")
    fig.show()      

    x = [110, 100, 130, 120, 150, 140, 170, 160, 180, 10, 0, 30, 20, 50, 40, 70, 60, 90, 80]
    y1 =[1.37041,0.932218,0.638898,0.319449,0.0310855,0.0001792,0.0413884,0.234315,0.605087,0.882136,1.02807,1.17223,1.13955,1.10105,0.990579,0.86636,0.761615,0.606906,0.477322]
    y2= [1.89103,2.2145,2.68132,3.13768,2.92068,2.53678,2.09305,1.35113,0.843754,0.391628,0.165586,0.0411888,0.000274262,0.0145772,0.063234,0.136357,0.217084,0.305137,0.389801]
    y1 = [1.6703567225501708, 1.386095070448514, 1.43737947280643, 1.34881376139331, 0.85723783413488663, 0.9158121288519564, 0.48003806599699189, 0.73330085260910549, 0.3629509716814075, 1.5638544731578126, 1.9124178429348806, 0.27919070777954863, 0.79633529993298613, 0.050264821456170491, 0.015284645164238357, 0.73301110653992751, 0.30115149916910217, 1.1421106720870242, 1.11201948449689]
    y2 = [0.029124547114672415, 0.11948762506787053, 0.020782968186636077, 7.3023485871301033e-05, 0.11550297386198953, 0.056855085915262539, 0.24058634006911733, 0.19190448733391069, 0.30947459589708609, 1.9967531880167566, 1.4870048197596841, 2.2913084908629235, 2.1376457565912736, 2.1579134837902068, 1.3551588961339354, 1.1534692696104769, 1.5281755585580536, 0.29168951289039557, 0.66353410469396068]    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('log')
    ax.plot((z_angles*360/N),chi_2)
    ax.plot(x,y1, 'ro', x,y2,'go')
    ax.legend(['Eigen Chi1','Eigen Chi2','HFSS Chi1','HFSS Chi2'])
    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("Cross-Kerr [MHz]")
    ax.set_xlim([0,180])
    ax.set_ylim([.005,50])
    ax.set_title("Cross-kerr")
    fig.show()  
    
    anharm1 = [271.4736003897878, 266.18367918217245, 271.2273441134771, 266.63064998656864, 264.92007310057079, 261.79795937563932, 261.47769875530878, 266.19128202250874, 261.97889019954016, 267.03839226738683, 264.04993257807695, 264.7051767467965, 264.46293013397133, 267.94722275527516, 244.12017327432056, 266.72164876946766, 263.82055841093455, 262.14637849308991, 266.51096774286486]    
    anharm2 = [0.0025693949398417314, 0.0018044490464484745, 0.0019043615934063831, 0.0017058227955184591, 0.00069347019996602883, 0.00080091519558879681, 0.00022032141354985606, 0.00050502230609469326, 0.00012570994531672103, 0.0022895966310781139, 0.0034627371140250746, 7.3617233584560626e-05, 0.00059946956422031634, 2.3573226940369144e-06, 2.3924730867505281e-07, 0.00050361986436962543, 8.5941203746651233e-05, 0.0012439775010371333, 0.0011599779029486444]
    anharm3 = [7.8114339977888503e-07, 1.340924863258504e-05, 3.981270473101838e-07, 4.9998279352586791e-12, 1.2589586752358754e-05, 3.0868277221689016e-06, 5.5341036064818298e-05, 3.4587282478861394e-05, 9.1395269894376408e-05, 0.0037326311583905647, 0.0020935276449414592, 0.0049584358953834321, 0.0043196312790963285, 0.0043446901181138572, 0.0018806881147295641, 0.0012470785199420529, 0.002212981952430576, 8.1140518151845241e-05, 0.00041300130330546524]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('log')
    ax.plot((z_angles*360/N),alpha)
    ax.plot(x,anharm1, 'ro', x,anharm2,'go',x,anharm3, 'bo')
    ax.set_xlabel("Zmon Angle [deg]")
    ax.set_ylabel("Anharmonicity [MHz]")
    ax.set_xlim([0,180])
    ax.set_ylim([5E-8,5000])
    ax.set_title("Anharmonicity")
    fig.show()  

    print sum(np.abs(diff))/len(diff)
    #%%
def collect_and_perturb():
    N =100
    #One_Trial(100)
    R = 11
    z_angles = np.linspace(0,50,R)
    freq_eigen= []
    freq_hfss = []
    freq_eigen2= []
    diff = []
    for z_angle in z_angles:       
        design = hfss.get_active_design()
        design.set_variable('zmonAngle',(u'%.4f' % (z_angle*360/N)))
        f,m,d,f2 = One_Trial(N, collect_new_data = False, qubit = True, qu_theta = z_angle,verbose = False)#(360-z_angle)/360*N)
        print f,m,d,f2 
        freq_eigen.append(f)
        freq_eigen.append(f2)
        freq_hfss.append(m)
        diff.append(d)
    
    print "Z_angles", z_angles
    print "freq eigen", freq_eigen
    print freq_eigen2
    print "freq hfss", freq_hfss
    print "diff", diff
    
    plt.plot((z_angles*360/N),freq_eigen)
    plt.xlabel("Zmon Angle [deg]")
    plt.ylabel("Eigenmode Frequencies [GHz]")
    plt.title("Perturbatively Calculated Eigenmodes")
    plt.show()    
    plt.plot((z_angles*360/N),freq_eigen2)
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
    #N=100    
    #One_Trial(N, collect_new_data = False, qubit = True, qu_theta = 0*z_angle,verbose = True)
    #collect_and_perturb()
    multi_perturb()