# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 20:49:48 2016

@author: Devin
"""
import hfss, numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class simulated_wg(object):
    def __init__(self):
        prefix = "../data/parameters"
        name = prefix + "capacitance" + ".npy"
        self.capacitance = np.load(name)
        name = prefix + "voltage" + ".npy"
        self.voltage = np.load(name)
        name = prefix + "inductance" + ".npy"
        self.inductance = np.load(name)
        name = prefix + "current" + ".npy"
        self.current = np.load(name)
        name = "../data/parametersangles.npy"
        self.angles = np.load(name)
        name = "../data/parameterseigenmodes.npy"
        self.eigenmodes = np.load(name)
    
    def build_L_mat(self):
        pass

    def build_C_mat(self):
        pass

def interpolate_outliers(angle, data, m=1.5,mv_av = 8):
    '''
    Function to smooth outliers from the data set. Applys moving
    average smoothing and cyclic boundary conditions. Threshold
    is set by:
    m - number of standard deviations from average which defines outliers
    mv_av - number of points in each direction used for average
    '''
    data2 = []
    mean = np.mean(data)
    std = np.std(data)
    for i in range(len(data)):
        if abs(data[i] - mean) < m*std:
            data2.append(data[i])
        else:
            ang_last = len(angle)
            if (i-mv_av < 0):
                comb = np.concatenate((data[0: i+8],data[ang_last + i - mv_av: ang_last]))
                print comb
                data2.append(np.mean(comb))
            elif(i+mv_av > ang_last):
                comb = np.concatenate((data[i-8: ang_last], data[0:i+mv_av - ang_last]))
                print comb
                data2.append(np.mean(comb))
            else:
                data2.append(np.mean(data[i-8:i+8]))
    return data2


def main():
    wg = simulated_wg()
    voltage = interpolate_outliers(wg.angles,wg.inductance)
    plt.plot(wg.angles,voltage)
    plt.show()

    
if __name__ == "__main__":
    main()            
        
        