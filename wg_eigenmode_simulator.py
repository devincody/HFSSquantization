# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 20:49:48 2016
@author: Devin

Eigenmode simulator
Takes in LCVI arrays and builds LC arrays, solves for eigenmodes
"""
import sys
if ~(r'F:\Documents\Yale\Devoret_Research\HFSSpython\pyHFSS' in sys.path):
    sys.path.append(r'F:\Documents\Yale\Devoret_Research\HFSSpython\pyHFSS')
import hfss, numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
import math

class simulated_wg(object):
    def __init__(self, name_variable = 0):
        prefix = "../data/parameters"
        name = prefix + "capacitance" + str(name_variable) + ".npy"
        self.capacitance = np.load(name)
        name = prefix + "voltage" + str(name_variable) + ".npy"
        self.voltage = np.load(name)
        name = prefix + "inductance" + str(name_variable) + ".npy"
        self.inductance = np.load(name)
        name = prefix + "current" + str(name_variable) + ".npy"
        self.current = np.load(name)
        name = "../data/parametersangles" + str(name_variable) + ".npy"
        self.angles = np.load(name)
        name = "../data/parameterseigenmodes" + str(name_variable) + ".npy"
        self.eigenmodes = np.load(name)
    
    def build_L_mat(self, N = 100, qubit = False, qu_theta = 0, verbose = False):
        #smooth the data
        smooth_l = interpolate_outliers(self.angles, self.inductance, plot_me = verbose)  
        smooth_l = moving_average(smooth_l, int(6*len(self.angles)/100.), plot_me = verbose) 
        
        if N == -1:
            size = len(smooth_l)
        else:
            size = N
        d_l = 5.89E-3*2*np.pi/N ##average circumfrence divided by numb nodes
        
        if qubit == True:
            design = hfss.get_active_design()
            Ravg = float(design.get_variable_value('RavgR').split('mm')[0])*1E-3
            Rx = float(design.get_variable_value('Rxoff').split('mm')[0])*1E-3
            Rt = float(design.get_variable_value('Rthick').split('mm')[0])*1E-3
            d = float(design.get_variable_value('Wspacing').split('mil')[0])*2.54E-5 #meters
            
            Zperplen = float(design.get_variable_value('zPerpLen').split('um')[0])*1E-6
            Zperpwid = float(design.get_variable_value('zPerpWid').split('um')[0])*1E-6
            gaplen = float(design.get_variable_value('GapLength').split('um')[0])*1E-6
            gapwid = float(design.get_variable_value('GapWidth').split('um')[0])*1E-6            
            mu = 4*np.pi*1E-7
            
            print (qu_theta*360/N)
            Length = Ravg + Rt/2 - (Rx*np.cos(qu_theta*2*np.pi/N)+np.sqrt((Ravg - Rt/2)**2 - (Rx**2)*(np.sin(qu_theta*2*np.pi/N))**2))
            
            Length_reduced = Length - (Zperplen + 2*gaplen)
            Width = 2*np.pi*Ravg/N
            Width_reduced = Width - (Zperpwid + 2*gapwid)
            
            Ljj = 9E-9
            Ljj_side = mu*d*(Width_reduced/Length + (Zperpwid + 2*gapwid)/(Length_reduced))
            print Ljj_side            
            #Ljj_side = mu*d*(Width/Length_reduced)       
            #print Ljj_side  
            if verbose:
                print "Ljj_side", Ljj_side

            self.L = np.zeros((size+1,size+1))
            for i in range(size+1):
                if i < qu_theta:
                    #print (d_l*smooth_l[i])**-1
                    self.L[i][i]          +=  (d_l*smooth_l[i])**-1
                    self.L[(i-1)%(size+1)][i] += -(d_l*smooth_l[i])**-1
                    self.L[i][(i-1)%(size+1)] += -(d_l*smooth_l[i])**-1
                    self.L[(i-1)%(size+1)][(i-1)%(size+1)] +=  (d_l*smooth_l[i])**-1
                elif i > qu_theta+1:
                    self.L[i][i]          +=  (d_l*smooth_l[i-1])**-1
                    self.L[(i-1)%(size+1)][i] += -(d_l*smooth_l[i-1])**-1
                    self.L[i][(i-1)%(size+1)] += -(d_l*smooth_l[i-1])**-1
                    self.L[(i-1)%(size+1)][(i-1)%(size+1)] +=  (d_l*smooth_l[i-1])**-1
                elif i == qu_theta:
                    self.L[(i+1)%(size+1)][(i+1)%(size+1)] +=  Ljj_side**-1 #inductance of transmission line around island
                    self.L[(i-1)%(size+1)][(i+1)%(size+1)] += -Ljj_side**-1
                    self.L[(i+1)%(size+1)][(i-1)%(size+1)] += -Ljj_side**-1 
                    self.L[(i-1)%(size+1)][(i-1)%(size+1)] +=  Ljj_side**-1
                    
                    ##draw bridge
                    self.L[i][i]          +=  Ljj**-1
                    self.L[(i+1)%(size+1)][i] += -Ljj**-1
                    self.L[i][(i+1)%(size+1)] += -Ljj**-1
                    self.L[(i+1)%(size+1)][(i+1)%(size+1)] +=  Ljj**-1
                    
        else:
            self.L = np.zeros((size,size))
            for i in range(size):
                self.L[i][i]          +=  (d_l*smooth_l[i])**-1
                self.L[(i-1)%size][i] += -(d_l*smooth_l[i])**-1
                self.L[i][(i-1)%size] += -(d_l*smooth_l[i])**-1
                self.L[(i-1)%size][(i-1)%size] +=  (d_l*smooth_l[i])**-1
            
        if verbose:
            print "L:", self.L

    def build_C_mat_parallel(self, N =100, qubit = False, qu_theta = 0, verbose = False):
        #smooth the data       
        smooth_c = interpolate_outliers(self.angles, self.capacitance, plot_me = verbose)
        smooth_c = moving_average(smooth_c, int(6*len(self.angles)/100.), plot_me = verbose)        
        
        if N == -1:
            size = len(smooth_l)
        else:
            size = N
            
        d_l = 5.89*.001*2*np.pi/len(self.angles)  ##average circumfrence divided by numb nodes

        if qubit == True:
            design = hfss.get_active_design()
            Ravg = float(design.get_variable_value('RavgR').split('mm')[0])*1E-3
            Rx = float(design.get_variable_value('Rxoff').split('mm')[0])*1E-3
            Rt = float(design.get_variable_value('Rthick').split('mm')[0])*1E-3
            d = float(design.get_variable_value('Wspacing').split('mil')[0])*2.54E-5
            
            Zperplen = float(design.get_variable_value('zPerpLen').split('um')[0])*1E-6
            Zperpwid = float(design.get_variable_value('zPerpWid').split('um')[0])*1E-6
            gaplen = float(design.get_variable_value('GapLength').split('um')[0])*1E-6
            gapwid = float(design.get_variable_value('GapWidth').split('um')[0])*1E-6
            
            #define demsions of node in which qubit is embedded
            Length = Ravg + Rt/2 - (Rx*np.cos(qu_theta*2*np.pi/N)+np.sqrt((Ravg - Rt/2)**2 - Rx**2*np.sin(qu_theta*2*np.pi/N)**2))
            Length_reduced = Length - (Zperplen + 2*gaplen) #reduced = ring dimesion minus island dimension
            Width = 2*np.pi*Ravg/N #circumference divided by number of nodes
            Width_reduced = Width - (Zperpwid + 2*gapwid)
            epsilon = 8.85E-12            
            
            self.C = np.zeros((size+1,size+1))
            Cjj_side = epsilon*(Length*Width_reduced/2+Length_reduced*(Zperpwid + 2*gapwid)/2)/d #
            Cjj = epsilon*(Zperplen*Zperpwid)/d #capacitance of island
            Cjj_spacing = 1*6.487E-14 #capacitance between a shore and island
            if verbose:            
                print "Cjj_side", Cjj_side, "Cjj_spacing", Cjj_spacing
            for i in range(size+1):
                if i < qu_theta:
                    self.C[i][i] += (d_l*smooth_c[i]/2)
                    self.C[(i-1)%(size+1)][(i-1)%(size+1)] += (d_l*smooth_c[i]/2)
                elif i > qu_theta+1:
                    self.C[i][i] += (d_l*smooth_c[i-1]/2)
                    self.C[(i-1)%(size+1)][(i-1)%(size+1)] += (d_l*smooth_c[i-1]/2)
                elif i == qu_theta:
                    self.C[i][i] += Cjj_spacing + Cjj # Qu_theta denotes where the qubit is
                    self.C[(i-1)%(size+1)][i] += -Cjj_spacing/2 #couple qubit to previous node
                    self.C[i][(i-1)%(size+1)] += -Cjj_spacing/2
                    self.C[(i-1)%(size+1)][(i-1)%(size+1)] +=  Cjj_spacing/2
                    
                    self.C[(i+1)%(size+1)][i] += -Cjj_spacing/2#couple qubit to next node
                    self.C[i][(i+1)%(size+1)] += -Cjj_spacing/2
                    self.C[(i+1)%(size+1)][(i+1)%(size+1)] +=  Cjj_spacing/2
                elif i == qu_theta+1:
                    self.C[(i)%(size+1)][(i)%(size+1)] +=  Cjj_side
                    self.C[(i-2)%(size+1)][(i-2)%(size+1)] += Cjj_side
        else:
            self.C = np.zeros((size,size))
            for i in range(size):
                self.C[i][i] += (d_l*smooth_c[i]/2)
                self.C[(i-1)%size][(i-1)%size] += (d_l*smooth_c[i]/2)
                
        if verbose:
            print "C:", self.C

    def build_C_mat(self, verbose = False):
        #depreciated
        smooth_c = interpolate_outliers(self.angles, self.capacitance, plot_me = verbose)
        smooth_c = moving_average(smooth_c, int(6*len(self.angles)/100.), plot_me = verbose)        
        size = len(smooth_c)
        self.C = np.zeros((size,size))
        d_l = 5.89*.001*2*np.pi/len(self.angles)  ##average circumfrence divided by numb nodes
        for i in range(size):
            self.C[i][i] = (d_l*smooth_c[i])
        if verbose:
            print "C:", self.C
            
    def build_L_mat_test(self,verbose = False):
        #depreciated
        smooth_l = interpolate_outliers(self.angles, self.inductance)   
        smooth_l = interpolate_outliers(self.angles, smooth_l)
        size = len(smooth_l)
        self.L = np.zeros((size,size))
        d_l = 0.001#5.89*.001*2*np.pi/100
        #v = 1;
        for i in range(size):
            v                     = (d_l*(.4+.05*.4*np.random.randn())*np.pi*10**(-7))**-1
            self.L[i][i]          +=  v
            self.L[(i+1)%size][i] += -v
            self.L[i][(i+1)%size] += -v
            self.L[(i+1)%size][(i+1)%size] +=  v
        print "L:", self.L

    def build_C_mat_test(self):
        #depreciated
        smooth_c = interpolate_outliers(self.angles, self.capacitance)
        size = len(smooth_c)
        self.C = np.zeros((size,size))
        d_l = 0.001#5.89*.001*2*np.pi/100
        for i in range(size):
            self.C[i][i] = (d_l*(10+.05*10*np.random.randn())*8.8541878176*10**(-12))
            #self.C[(i+1)%size][i] += -(d_l*2.3*10**(-10))
            #self.C[i][(i+1)%size] += -(d_l*2.3*10**(-10))
            #self.C[(i+1)%size][(i+1)%size] += (d_l*2.3*10**(-10))
        print "C:", self.C

    def test_interpolate(self,plot_me = True):
        potted = self.capacitance
        voltage = interpolate_outliers(self.angles,potted, plot_me = False)
        voltage = interpolate_outliers(self.angles,voltage, plot_me = False)
        if plot_me:
            plt.plot(self.angles,voltage, self.angles, potted)
            ax = plt.gca()
            ax.set_ylim(min(voltage)-abs(min(voltage)*.05), max(voltage)*1.05)
            plt.title("Test of Interpolation")            
            plt.show()
        
    def get_frequencies(self,qu_theta=0, verbose = False):
        n_modes = 3
        LC = np.dot(np.linalg.inv(self.C),self.L)
        w,v = np.linalg.eig(LC)
        idx = w.argsort()[::1]
        w = w[idx]
        w = np.sqrt(w)/(2*np.pi) #obtain eigen frequencies, not angular velocities
        v = v[:,idx] #sort eigenvectors
        self.calc_freq = w

        if verbose:
            plt.figure()
            plt.plot(np.sort(self.calc_freq)[0:8]/(10**9))
            plt.title("Frequency Eigenmodes")
            plt.xlabel("Eigenmode index")
            plt.ylabel("Frequency (Hz)")
            plt.show()
            print "inv C"
            print np.linalg.inv(self.C)
            print "LC"
            print LC
            print "frequencies:"
            print self.calc_freq
            print "eigenvalues: ", w
            print "eigenvectors: ", v      
        
        Ljj = 9E-9 #inductance of josephson junction
        h = 6.626E-34 #plank's constant
        hbar= h/(2*np.pi)
        e = 1.602E-19 #electron charge in coulombs
        reduced_flux_quant = hbar/(2*e)
        EJ = reduced_flux_quant**2/Ljj #Josephson energy in junction 
        
        i = 0
        if self.calc_freq[0]< 4E9 or math.isnan(self.calc_freq[0]):
            i = 1
            print "first mode is garbage"
            
#        E1 = .5*np.dot(v[i,:].T,np.dot(self.L,v[i,:]))
#        Ejj = (v[i,qu_theta+1]-v[i,qu_theta])**2/Ljj
        E = np.zeros(n_modes)
        Ejj = np.zeros(n_modes)
        EPR = np.zeros(n_modes)
        nu = np.zeros(n_modes)
        chi = np.zeros((n_modes,n_modes))
            
        for r in range(3):
            E[r] = .5*np.dot(v[:,r+i].T,np.dot(self.L,v[:,r+i])) # Inductive energy in mode r
            Ejj[r] = .5*(v[qu_theta,r+i]-v[qu_theta-1,r+i])**2/Ljj # Energy stored in junction inductance
            EPR[r] = Ejj[r]/E[r] # Definition of Energy participation ratio
            nu[r] = w[r+i] # frequency of rth node
        
        if verbose:
            print "EPR: ", EPR
            
        for m in range(n_modes):
            for n in range(n_modes):
                chi[n][m] = nu[n]*nu[m]*EPR[m]*EPR[n]*h*np.pi/(2*EJ)
        if verbose:
            print "CHI: "
            print chi/1E6 #report in MHz
            print "frequencies: ",nu
            print "Joshpson Energy: ",EJ

        return np.sort(self.calc_freq),chi
        

def interpolate_outliers(angle, data, threshold=.5, window = 12, plot_me = False):
    '''
    Function to smooth outliers from the data set. Applys moving
    average smoothing and cyclic boundary conditions. Threshold
    is set by:
    threshold - number of standard deviations from average which defines outliers
    window - number of points in each direction used for average
    '''
    df = pd.DataFrame({'parameter':data},index=angle)
    #mean_data = np.mean(df['parameter'])
    df['data_mean'] = pd.rolling_median(df['parameter'].copy(), window=window, center=True).fillna(method='bfill').fillna(method='ffill')
    difference = np.abs(df['parameter'] - df['data_mean'])#mean_data)#
    outlier_idx = difference > threshold*df['parameter'].std()    
    #df['data_mean'].plot()
#    s = df['parameter'].copy()
#    s[outlier_idx] = np.nan
#    s.interpolate(method='spline', order=1, inplace=True)
#    df['cleaned_parameter'] = s
    tst = np.array(outlier_idx)
    datamean= np.array(df['data_mean'])
    s = np.array(df['parameter'])
    itms = len(outlier_idx)
    for i in range(itms):
        if tst[i] == True or tst[(i-1)%itms] == True or tst[(i+1)%itms] == True or tst[(i-2)%itms] == True or tst[(i+2)%itms] == True:
            tmp =  datamean[i]
            s[i] = tmp
    #print s
    df['cleaned_parameter'] = s
    
    if (plot_me == True):
        figsize = (7, 2.75)
        fig, ax = plt.subplots(figsize=figsize)
        df['parameter'].plot(title = "cleaned vs unclean Parameter")
        df['cleaned_parameter'].plot()
        ax.set_ylim(min(df['cleaned_parameter']), max(df['cleaned_parameter']))
    return np.array(df['cleaned_parameter'])

def moving_average(x, n, plot_me):
    y = np.zeros(len(x))
    for i in range(len(x)):
        summ = 0
        for j in range(n):
            summ += x[int(i-n/2+j)%len(x)]
        y[i] = summ/n
    if (plot_me == True):
        figsize = (7, 2.75)
        fig, ax = plt.subplots(figsize=figsize)
        plt.plot(x, 'g', y,'r')
        plt.title("moving average")
        plt.show()
    return y
        
def main():
    sim_wg = simulated_wg()
    sim_wg.test_interpolate()
    sim_wg.build_L_mat_test()
    sim_wg.build_C_mat_test()
    sim_wg.get_frequencies(verbose = True)
    
if __name__ == "__main__":
    main()            