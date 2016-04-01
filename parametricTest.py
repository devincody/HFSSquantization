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

class waveguide(object):
    def __init__(self, angle_s = 0, angle_e = 360, angle_n = 100):
        self.proj = get_active_project()
        self.design = get_active_design()
        self.name = self.design.get_setup_names()
        self.setup = self.design.get_setup(self.name[0])
        self.angles = np.linspace(angle_s, angle_e, angle_n)
        self.capacitance = []
        self.voltage = []
        self.inductance = []
        self.current = []
        try:
            self.solutions = self.setup.get_solutions()
            self.eigenmodes = self.solutions.eigenmodes()
        except:
            print "Analyzing Geometry... this might take a while"
            self.setup.analyze()
            self.solutions = self.setup.get_solutions()
            self.eigenmodes = self.solutions.eigenmodes()
    
    def calc_voltage(self, fields, line = "intLineVolt"):
        '''Function to calculate WaveGuide Voltage
        line = integration line between plates'''
        self.design.Clear_Field_Clac_Stack()
        comp = fields.Vector_E
        exp = comp.integrate_line_tangent(line)
        V = exp.evaluate(phase = 0) 
        self.design.Clear_Field_Clac_Stack()
        return V
    
    def calc_current(self, fields, line = "IntLineCurrent"):
        '''Function to calculate WaveGuide Current
        line = integration line between plates'''
        self.design.Clear_Field_Clac_Stack()
        comp = fields.Vector_H
        exp = comp.integrate_line_tangent(line)
        I = exp.evaluate(phase = 90)
        self.design.Clear_Field_Clac_Stack()
        return I
    
    def calc_inductance(self, fields, surf = "CrossSecIntSurf", line = "IntLineCurrent"):
        '''Function to calculate WaveGuide Inductance
        surf = integration surface between plates    
        line = integration line between plates
        returns current as secondary parameter'''  
        self.design.Clear_Field_Clac_Stack()
        I = self.calc_current(fields, line)
        mu = 4*np.pi*10**-7
        Mag_H_Sq = fields.Mag_H ** 2
        Surf_H = Mag_H_Sq.integrate_surf(surf)
        preinductance = Surf_H.evaluate(phase = 90)
        L = preinductance*mu/(I**2)   
        self.design.Clear_Field_Clac_Stack()
        return L, I
        
    def calc_capacitance(self,fields, surf = "CrossSecIntSurf", line = "intLineVolt"):
        '''Function to calculate WaveGuide Voltage
        surf = integration surface between plates 
        line = integration line between plates
        returns voltage as secondary parameter''' 
        self.design.Clear_Field_Clac_Stack()
        V = self.calc_voltage(fields, line)   
        epsilon = 8.8541878176*10**(-12)    
        Mag_E_Sq = fields.Mag_E ** 2
        Surf_E = Mag_E_Sq.integrate_surf(surf)
        precapacitance = Surf_E.evaluate()
        C = precapacitance*epsilon/(V**2)
        self.design.Clear_Field_Clac_Stack()
        return C, V

    def compute_LCVI(self, verbose = False):
        print "Compute LCVI"
        for i in self.angles:
            self.design.set_variable('th',(u'%.2fdeg' % (i)))
            fields = self.setup.get_fields()
            C, V = self.calc_capacitance(fields)
            L, I = self.calc_inductance(fields)
            self.capacitance.append(C)
            self.voltage.append(V)
            self.inductance.append(L)
            self.current.append(I)
            if (verbose):
                print "#################"
                print "Angle: " , i
                print "voltage:", V
                print "current:", I
                print "capacitance:", C
                print "inductance:", L
    
    def plot(self,scale_factor = 1.2):
        labels = ["Capacitance", "Voltage", "Inductance", "Current"]
        ylabel_name = "Angle around disk (Degrees)"
        
        ang, cap = reject_outliers(self.angles, self.capacitance)
        plt.plot(ang, cap)
        plt.title(labels[0])
        plt.xlabel(ylabel_name)
        plt.ylabel(labels[0] + "(Farads per Meter)")
        plt.axis([0,360, min(cap), max(cap)])
        plt.show()
        
        ang, vot = reject_outliers(self.angles,self.voltage)
        plt.plot(ang, vot) 
        plt.title(labels[1])
        plt.xlabel(ylabel_name)
        plt.ylabel(labels[1] + "(Volts)")
        plt.axis([0,360, min(vot)*scale_factor, max(vot)*scale_factor])
        plt.show()
        
        ang, ind = reject_outliers(self.angles,self.inductance)
        plt.title(labels[2])
        plt.xlabel(ylabel_name)
        plt.ylabel(labels[2] + "(Henries per Meter)")
        plt.axis([0,360, min(ind)*scale_factor, max(ind)*scale_factor])
        plt.plot(ang, ind)
        plt.show()
        
        ang, cur = reject_outliers(self.angles,self.current)
        plt.title(labels[3])
        plt.xlabel(ylabel_name)
        plt.ylabel(labels[3] + "(Amps)")
        plt.axis([0,360, min(cur)*scale_factor, max(cur)*scale_factor])
        plt.plot(ang, cur)
        plt.show()
        
    def save(self, name_variable = 0):
        prefix = "../data/parameters"
        labels = ["capacitance", "voltage", "inductance", "current"]
        for i in labels:
            name = prefix + i + str(name_variable) + ".npy"
            ii = "self." + i
            np.save(name, eval(ii))
        name = "../data/parametersangles"+ str(name_variable) + ".npy"
        np.save(name, self.angles)
        name = "../data/parameterseigenmodes"+ str(name_variable) + ".npy"
        np.save(name, self.eigenmodes)
        print self.capacitance        
        
    def load(self, name_variable = 0):
        prefix = "../data/parameters"
        #labels = ["capacitance", "voltage", "inductance", "current"]
        name = prefix + "capacitance"+ str(name_variable) + ".npy"
        self.capacitance = np.load(name)
        name = prefix + "voltage"+ str(name_variable) + ".npy"
        self.voltage = np.load(name)
        name = prefix + "inductance"+ str(name_variable) + ".npy"
        self.inductance = np.load(name)
        name = prefix + "current"+ str(name_variable) + ".npy"
        self.current = np.load(name)
        
    def set_scalex(self, value):
        self.design.set_variable('scalex',(u'%.4f' % (value)))
        
    def set_scalez(self, value):
        self.design.set_variable('scalez',(u'%.4f' % (value)))

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

def main():
    wg = waveguide() 
    wg.compute_LCVI()      
    print wg.inductance
    wg.save()
    wg.load()   
    print wg.inductance
    wg.plot()
    hfss.release()
    
if __name__ == "__main__":
    main()    
