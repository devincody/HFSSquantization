# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:38:28 2016

@author: Devin Cody
"""
#from hfss import *
import sys
if ~(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS' in sys.path):
    sys.path.append(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS')
import hfss,numpy as np
import matplotlib.pyplot as plt
from hfss import get_active_design
from hfss import get_active_project

def calc_voltage(fields, line = "intLineVolt"):
    '''Function to calculate WaveGuide Voltage
    line = integration line between plates'''
    comp = fields.Vector_E
    exp = comp.integrate_line_tangent(line)
    voltage = exp.evaluate(phase= 0) 
    design.Clear_Field_Clac_Stack()
    return voltage

def calc_current(fields, line = "IntLineCurrent"):
    '''Function to calculate WaveGuide Current
    line = integration line between plates'''
    comp = fields.Vector_H
    exp = comp.integrate_line_tangent(line)
    current = exp.evaluate(phase= 90)
    design.Clear_Field_Clac_Stack()
    return current

def calc_inductance(fields, surf = "CrossSecIntSurf", line = "IntLineCurrent"):
    '''Function to calculate WaveGuide Inductance
    surf = integration surface between plates    
    line = integration line between plates
    returns current as secondary parameter'''   
    mu = 4*np.pi*10**-7
    i = calc_current(fields, line)
    Mag_H_Sq = fields.Mag_H ** 2
    Surf_H = Mag_H_Sq.integrate_surf(surf)
    preinductance = Surf_H.evaluate()
    inductance = preinductance*mu/(i**2)   
    design.Clear_Field_Clac_Stack()
    return inductance, i
    
def calc_capacitance(fields, surf = "CrossSecIntSurf", line = "intLineVolt"):
    '''Function to calculate WaveGuide Voltage
    surf = integration surface between plates 
    line = integration line between plates
    returns voltage as secondary parameter'''      
    epsilon = 8.8541878176*10**(-12)    
    v = calc_voltage(fields, line)    
    Mag_E_Sq = fields.Mag_E ** 2
    Surf_E = Mag_E_Sq.integrate_surf(surf)
    precapacitance = Surf_E.evaluate()
    capacitance = precapacitance*epsilon/(v**2)
    design.Clear_Field_Clac_Stack()
    return capacitance, v

def reject_outliers(angle, data, m=2):
    angle2 = []
    data2 = []
    mean = np.mean(data)
    std = np.std(data)   
    for i in range(len(data)):
        print i, data
        if abs(data[i] - mean) < m*std:
            data2.append(data[i])
            angle2.append(angle[i])
    return angle2, data2

def main():
    proj = get_active_project()
    design = get_active_design()
    name = design.get_setup_names()
    setup = design.get_setup(name[0])
    #setup.analyze()
    
    
    angles = np.linspace(0,360,100)
    capacitance = []
    voltage = []
    inductance = []
    current = []
    for i in angles:
        design.set_variable('th',(u'%.2fdeg' % (i)))
        fields = setup.get_fields()
        C,V = calc_capacitance(fields)
        L,I = calc_inductance(fields)
        capacitance.append(C)
        voltage.append(V)
        inductance.append(L)
        current.append(I)
        print "#################"
        print "Angle: " , i
        print "voltage:", V
        print "current:", I
        print "capacitance:", C
        print "inductance:", L

def plot():
    ang, cap = reject_outliers(angles,capacitance)
    plt.plot(ang, cap)
    plt.show()
    
    ang, vot = reject_outliers(angles,voltage)
    plt.plot(ang, vot) 
    plt.show()
    
    ang, ind=reject_outliers(angles,inductance)
    plt.plot(ang, ind)
    plt.show()
    
    ang, cur = reject_outliers(angles,current)
    plt.plot(ang, cur)
    plt.show()
    
def save():
    prefix = "../data/parameters"
    labels = ["capacitance", "voltage", "inductance", "current"]
    for i in labels:
        name = prefix + i + ".npy"
        np.save(name, eval(i))

def load():
    prefix = "../data/parameters"
    labels = ["capacitance", "voltage", "inductance", "current"]
    
    name = prefix + "capacitance" + ".npy"
    capacitance = np.load(name)
    
    name = prefix + "voltage" + ".npy"
    voltage = np.load(name)
    
    name = prefix + "inductance" + ".npy"
    inductance = np.load(name)
    
    name = prefix + "current" + ".npy"
    current = np.load(name)
    
if __name__ == "__main__":
    load()
    plot()
    hfss.release()