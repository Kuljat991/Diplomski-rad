# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 22:32:44 2018

@author: uljo
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import subprocess
import shutil
import math
from scipy.optimize import bisect
from scipy import interpolate
from math import pi, cos, sin

def Sigme_i_sve (X, drag, lift):
    h=10.
    def ogranicenje (qx, qy, x_max, y_max, Ix, Iy):
        Mx = qx * h**2 /2.
        My = qy * h**2 /2.
        sigma_x = Mx * y_max/Ix
        sigma_y = My * x_max/Iy
        return sigma_x , sigma_y
    
    def func_max_x(x,xt):
        delta_x=[]
        for i in range(len(x)):
            delta_x.append(abs(x[i]-xt))
        return max(delta_x)
    
    #posto je profil simetrican racuna se samo gornja krivulja od 0 do 180 rx ide s lijeva
    r_temp = X
    fi_temp = np.linspace(0., pi, len(r_temp))
    #print(len(fi_temp), len(r_temp))
    interp_r = interpolate.interp1d(fi_temp,r_temp)#,kind='cubic')
    
    n_points= 100
    fi_for_integrate=np.linspace(0.,pi,n_points)
    r_for_integrate=interp_r(fi_for_integrate)
    #plt.plot(fi_for_integrate,r_for_integrate)
    
    
    x_all=[]
    y_all=[]
    x_prvi_kvadr=[]
    y_prvi_kvadr=[]
    x_drugi_kvadr=[]
    y_drugi_kvadr=[]
    for i in range (len(fi_for_integrate)):
        if (fi_for_integrate[i]<=pi/2):
            x_prvi_kvadr.append(r_for_integrate[i]*cos(fi_for_integrate[i]))
            y_prvi_kvadr.append(r_for_integrate[i]*sin(fi_for_integrate[i]))
            x_all.append(x_prvi_kvadr[-1])
            y_all.append(y_prvi_kvadr[-1])
        else:
            x_drugi_kvadr.append(r_for_integrate[i]*cos(fi_for_integrate[i]))
            y_drugi_kvadr.append(r_for_integrate[i]*sin(fi_for_integrate[i]))
            x_all.append(x_drugi_kvadr[-1])
            y_all.append(y_drugi_kvadr[-1])
        plt.plot(x_prvi_kvadr,y_prvi_kvadr,'o')
        plt.plot(x_drugi_kvadr,y_drugi_kvadr,'o')           
    
    # Moment inercije Ix, Iy
    def r_phi(phi):
        phi=np.abs(phi)
        phi=phi%(2*pi)
        if(phi<pi):
            return interp_r(phi)
        else:
            return interp_r(2*pi-phi)
    
    ############################
    def integrand_x(fi_kut):
        return 1/4.*(r_phi(fi_kut))**4*(sin(fi_kut))**2
    fi_fi=np.linspace(0.,2*pi,101)    
    r_x=[]
    for i in range(len(fi_fi)):
        r_x.append(integrand_x(fi_fi[i]))
    Ix = integrate.simps(r_x,fi_fi)
    print(Ix)

    R=[]
    for i in range (len(fi_fi)):
        R.append((r_phi(fi_fi[i]))**2/2)
    A=integrate.simps(R,fi_fi)
    print(A)
    
    def integrand_y(fi_kut):
        return 1/4.*(r_phi(fi_kut))**4*(cos(fi_kut))**2
    r_y=[]
    for i in range(len(fi_fi)):
        r_y.append(integrand_y(fi_fi[i]))
    Iy = integrate.simps(r_y,fi_fi)
    print(Iy)
    
    def integrand_Ty(fi_kut):
        return 1/2.*(r_phi(fi_kut))**2*sin(fi_kut)
    r_Ty=[]
    for i in range(len(fi_fi)):
        r_Ty.append(integrand_Ty(fi_fi[i]))
    Ty = integrate.simps(r_Ty,fi_fi)
    #print(Ty)
    
    def integrand_Tx(fi_kut):
        return 1/2.*(r_phi(fi_kut))**2*cos(fi_kut)
    r_Tx=[]
    for i in range(len(fi_fi)):
        r_Tx.append(integrand_Tx(fi_fi[i]))
    Tx = integrate.simps(r_Tx,fi_fi)
    #print(Tx)
    x_max = func_max_x(x_all,Tx)
    y_max = max(y_all)
    #print('Tx='+str(Tx)+' Ty='+str(Ty))
    #drag stvara moment oko y-osi
    #lift stvara moment oko x-osi
    
    M_x = lift * h**2 /2.
    M_y = drag * h**2 /2.
    sigma_x=[]
    sigma_y=[]
    for i in range (len(x_prvi_kvadr)):
        sigma_x.append(- M_x / Ix * abs(y_prvi_kvadr[i]))
        sigma_y.append(- M_y / Iy * abs(x_prvi_kvadr[i]))
    for i in range (len(x_drugi_kvadr)):
        sigma_x.append(- M_x / Ix * abs(y_drugi_kvadr[i]))
        sigma_y.append( M_y / Iy * abs(x_drugi_kvadr[i]))  
    
    #fi_for_plot=np.rad2deg(fi_for_integrate)
    
    sigma_ekv_plus=[]
    sigma_ekv_minus=[]
    for i in range(len(sigma_x)):
        sigma_ekv_plus.append(sigma_x[i]+sigma_y[i])
        sigma_ekv_minus.append(-sigma_x[i]+sigma_y[i])
    
    sigma_ekv_plus_max=max([abs(max(sigma_ekv_plus)),abs(min(sigma_ekv_plus))])
    sigma_ekv_minus_max=max([abs(max(sigma_ekv_minus)),abs(min(sigma_ekv_minus))])
#    print(sigma_ekv_plus_max)
#    print(sigma_ekv_minus_max)
    
#    plt.plot(fi_for_plot, sigma_x, label='sigma_lift')
#    plt.plot(fi_for_plot, sigma_y, label='sigma_drag')
#    plt.plot(fi_for_plot, sigma_ekv_plus, label='sigma_superpon_plus')
#    plt.plot(fi_for_plot, sigma_ekv_minus, label='sigma_superpon_minus')
#    plt.xlabel('fi')
#    plt.ylabel('Sigma')
#    plt.legend()
#    plt.savefig('sigme.png')
    
    #Sigma_x, Sigma_y = ogranicenje (drag, lift, x_max, y_max, Ix, Iy)
    return 5, 5, Area, Ix, Iy
    
if __name__ == "__main__":
    a=12
    b=12
    fi_diag=math.atan2(b/2,a/2)
    #print(np.rad2deg(fi_diag))
    fi=np.linspace(0,pi,3610)
    r=[]
    for i in range (len(fi)):
        if (fi[i]<=fi_diag):
            x=a/2
            r.append(x/np.cos(fi[i]))
        elif (fi[i]>fi_diag and fi[i]<=pi-fi_diag):
            y=b/2
            r.append(y/np.sin(fi[i]))
        else:
            x=a/2
            r.append(abs(x/np.cos(fi[i])))
    #plt.plot(fi,r)
    
#    FileName = './proracun/rezultati_plot.txt'
#    text_X=np.loadtxt(FileName, skiprows=2, delimiter=',')
#    FileForces='./proracun/best/Best_iter/Forces_max.txt'
#    Forces=np.loadtxt(FileForces, skiprows=1, delimiter=',')
#    drag_all=Forces[:,1]
#    lift_all=Forces[:,2]
    #X=[2,2,2,2,2,2,2,2]
    X=r
    sigma_x, sigma_y, A, Ix, Iy = Sigme_i_sve(X,3000,3000)
    #print(sigma_x, sigma_y, A, Ix, Iy)
    
#    svasta=[]
#    for i in range(len(text_X)):
#        svasta.append(Sigme_i_sve(text_X[i],drag_all[i],lift_all[i]))
#
#    rezultati = open('./proracun/best/Best_iter/Sve_iter.txt', 'r')
#    contents = rezultati.readlines()
#    rezultati.close()
#    for i in range (len(svasta)):
#        contents.append(str(svasta[i][0])+'\t'+str(svasta[i][1])+'\t'+str(svasta[i][2])+'\t'+str(svasta[i][3])+'\t'+str(svasta[i][4])+'\t'+'\n')
#    out = open('./proracun/best/Best_iter/Sve_iter.txt', 'w')
#    for i in range(len(contents)):
#        out.writelines(str(contents[i]))
#    out.close()
