# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:38:28 2016

@author: Devin Cody
"""
#from hfss import *
import sys
if ~(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS' in sys.path):
    sys.path.append(r'F:\Documents\Yale\Junior Year\HFSSpython\pyHFSS')
import hfss, numpy as np
import matplotlib.pyplot as plt
from hfss import get_active_design
from hfss import get_active_project

def calc_voltage(design, fields, line = "intLineVolt"):
    '''Function to calculate WaveGuide Voltage
    line = integration line between plates'''
    comp = fields.Vector_E
    exp = comp.integrate_line_tangent(line)
    voltage = exp.evaluate(phase= 0) 
    design.Clear_Field_Clac_Stack()
    return voltage

def calc_current(design, fields, line = "IntLineCurrent"):
    '''Function to calculate WaveGuide Current
    line = integration line between plates'''
    comp = fields.Vector_H
    exp = comp.integrate_line_tangent(line)
    current = exp.evaluate(phase= 90)
    design.Clear_Field_Clac_Stack()
    return current

def calc_inductance(design, fields, surf = "CrossSecIntSurf", line = "IntLineCurrent"):
    '''Function to calculate WaveGuide Inductance
    surf = integration surface between plates    
    line = integration line between plates
    returns current as secondary parameter'''   
    mu = 4*np.pi*10**-7
    i = calc_current(design, fields, line)
    Mag_H_Sq = fields.Mag_H ** 2
    Surf_H = Mag_H_Sq.integrate_surf(surf)
    preinductance = Surf_H.evaluate()
    inductance = preinductance*mu/(i**2)   
    design.Clear_Field_Clac_Stack()
    return inductance, i
    
def calc_capacitance(design, fields, surf = "CrossSecIntSurf", line = "intLineVolt"):
    '''Function to calculate WaveGuide Voltage
    surf = integration surface between plates 
    line = integration line between plates
    returns voltage as secondary parameter'''      
    epsilon = 8.8541878176*10**(-12)    
    v = calc_voltage(design, fields, line)    
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
        if abs(data[i] - mean) < m*std:
            data2.append(data[i])
            angle2.append(angle[i])
    return angle2, data2

def init():
    proj = get_active_project()
    design = get_active_design()
    name = design.get_setup_names()
    setup = design.get_setup(name[0])
    #setup.analyze()
    return proj, design, setup

def analyze(proj, design, setup):
    angles = np.linspace(0,360,100)
    capacitance = []
    voltage = []
    inductance = []
    current = []
    for i in angles:
        design.set_variable('th',(u'%.2fdeg' % (i)))
        fields = setup.get_fields()
        C,V = calc_capacitance(design, fields)
        L,I = calc_inductance(design, fields)
        capacitance.append(C)
        print "capacitance length", len(capacitance)
        voltage.append(V)
        inductance.append(L)
        current.append(I)
        print "#################"
        print "Angle: " , i
        print "voltage:", V
        print "current:", I
        print "capacitance:", C
        print "inductance:", L
        
    return capacitance, voltage, inductance, current, angles

def plot(capacitance, voltage, inductance, current, angles):
    labels = ["Capacitance", "Voltage", "Inductance", "Current"]
    ylabel_name = "Angle around disk (Degrees)"
    
    ang, cap = reject_outliers(angles,capacitance)
    plt.plot(np.array(ang), np.transpose(cap))
    plt.title(labels[0])
    plt.xlabel(ylabel_name)
    plt.ylabel(labels[0] + "(Farads per Meter)")
    plt.axis([0,360, min(cap), max(cap)])
    plt.show()
    
    ang, vot = reject_outliers(angles,voltage)
    plt.plot(ang, vot) 
    plt.title(labels[1])
    plt.xlabel(ylabel_name)
    plt.ylabel(labels[1] + "(Volts)")
    plt.axis([0,360, min(vot), max(vot)])
    plt.show()
    
    ang, ind = reject_outliers(angles,inductance)
    plt.title(labels[2])
    plt.xlabel(ylabel_name)
    plt.ylabel(labels[2] + "(Henries per Meter)")
    plt.axis([0,360, min(ind), max(ind)])
    plt.plot(ang, ind)
    plt.show()   
    
    ang, cur = reject_outliers(angles,current)
    plt.title(labels[3])
    plt.xlabel(ylabel_name)
    plt.ylabel(labels[3] + "(Amps)")
    plt.axis([0,360, min(cur), max(cur)])
    plt.plot(ang, cur)
    plt.show()
    
def save(capacitance, voltage, inductance, current):
    prefix = "../data/parameters"
    labels = ["capacitance", "voltage", "inductance", "current"]
    for i in labels:
        name = prefix + i + ".npy"
        np.save(name, eval(i))
    print "save", len(capacitance)

def load():
    prefix = "../data/parameters"
    labels = ["capacitance", "voltage", "inductance", "current"]
    
    name = prefix + "capacitance" + ".npy"
    capacitance = np.load(name)
    print "load", len(capacitance)
     
    name = prefix + "voltage" + ".npy"
    voltage = np.load(name)
    
    name = prefix + "inductance" + ".npy"
    inductance = np.load(name)
    
    name = prefix + "current" + ".npy"
    current = np.load(name)
    return capacitance, voltage, inductance, current
    
def main():
    proj, design, setup = init()
    #capacitance, voltage, inductance, current, angles= analyze(proj, design, setup)    
    #save(capacitance, voltage, inductance, current)
    capacitance, voltage, inductance, current = load()
    plot(capacitance, voltage, inductance, current, angles)
    hfss.release()
    
if __name__ == "__main__":
    main()    
