# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:55:25 2017

@author: uljo
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import subprocess
import shutil

#citanje sila
def line2dict(line):
    tokens_unprocessed = line.split()
    tokens = [x.replace(")","").replace("(","") for x in tokens_unprocessed]
    floats = [float(x) for x in tokens]
    data_dict = {}
    data_dict['time'] = floats[0]
    force_dict = {}
    force_dict['pressure'] = floats[1:4]
    force_dict['viscous'] = floats[4:7]
    force_dict['porous'] = floats[7:10]
    moment_dict = {}
    moment_dict['pressure'] = floats[10:13]
    moment_dict['viscous'] = floats[13:16]
    moment_dict['porous'] = floats[16:19]
    data_dict['force'] = force_dict
    data_dict['moment'] = moment_dict
    return data_dict

newFileName='best_for_start'
list_file = os.listdir('./'+newFileName+'/postProcessing/forcesIncompressible')

time = []
drag = []
lift = []
moment = []
for i in range (len(list_file)):
    forces_file = './'+newFileName+'/postProcessing/forcesIncompressible/'+list_file[i]+'/forces.dat'
    if not os.path.isfile(forces_file):
        print ("Forces file not found at "+forces_file)
        print ("Be sure that the case has been run and you have the right directory!")
        print ("Exiting.")
        sys.exit()
    with open(forces_file,"r") as datafile:
        for line in datafile:
            if line[0] == "#":
                continue
            data_dict = line2dict(line)
            time += [data_dict['time']]
            drag += [data_dict['force']['pressure'][0] + data_dict['force']['viscous'][0]]
            lift += [data_dict['force']['pressure'][1] + data_dict['force']['viscous'][1]]
            moment += [data_dict['moment']['pressure'][2] + data_dict['moment']['viscous'][2]]
    datafile.close()

lift_za_plot=[]
drag_za_plot=[]
for i in range (len(lift)):
    lift_za_plot.append(lift[i]+1000)
    drag_za_plot.append(drag[i]+1000)
    
plt.plot(time[10000:],lift[10000:], label='Sila uzgona')
#plt.savefig('lift.png')
#plt.clf()
#plt.plot(time[10000:],drag[10000:])
#plt.savefig('drag.png')
plt.xlabel('t [s]')
plt.ylabel('F [N]')
#plt.legend()