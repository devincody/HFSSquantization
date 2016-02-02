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
    
    def build_L_mat(self):
        pass

    def build_C_mat(self):
        pass

def interpolate_outliers(angle, data, m=2):
    data2 = []
    mean = np.mean(data)
    std = np.std(data)
    f = interp1d(angle, data)
    for i in range(len(data)):
        if abs(data[i] - mean) < m*std:
            data2.append(data[i])
        else:
            data2.append(int(f(angle[i])))
    return data2


def main():
    wg = simulated_wg()    

    
if __name__ == "__main__":
    main()            
        
        