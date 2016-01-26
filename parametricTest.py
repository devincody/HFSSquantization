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

proj = get_active_project()
design = get_active_design()
name = design.get_setup_names()
setup = design.get_setup(name[0])
#setup.analyze()


angles = np.linspace(0,360,20)
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

plt.plot(angles, capacitance)
plt.plot(angles, voltage)
plt.plot(angles, inductance)
plt.plot(angles, current)
plt.show()



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

hfss.release()