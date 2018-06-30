# -*- coding: utf-8 -*-
"""
Created on Tue May 15 02:12:18 2018

@author: uljo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def plot_profile(R, filename):
    
    
    PHI = np.linspace(0., np.pi, len(R))
    interp_r = interpolate.interp1d(PHI, R, kind='cubic')
        
    def r_phi(phi):
        phi=np.abs(phi)
        phi=phi%(2*np.pi)
        if(phi < np.pi):
            return interp_r(phi)
        else:
            return interp_r(2*np.pi - phi)
        
    PHI = np.linspace(0, 2 * np.pi, 101)
    R = np.array([r_phi(phi) for phi in PHI])
    plt.figure(figsize=(10,6))
    plt.plot(R * np.cos(PHI), R * np.sin(PHI))
    plt.axis('equal')
    plt.xlim(-10, 10)
    plt.xlim(-7, 6)
    plt.xlabel('$x [m]$')
    plt.ylabel('$y [m]$')
    plt.savefig(filename)
    plt.close()
    
if __name__ == "__main__":
    R=[3.2207603893825074, 0.90801560199401887, 3.157813375761541, 2.7724328583924893, 2.0377827745169523, 1.8084522673720427, 1.9963162071726572, 6.0]
    #R=[3.69567485,1.11910148,4.08031717,2.60037542,1.05344528,2.62814699,5.04284431,1.58244596]
    profilepng = './profile_dvi.png'
    plot_profile(R, profilepng)