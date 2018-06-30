# -*- coding: utf-8 -*-
#solver = icoFoam
#foamMonitor -l -r 2 -i 5 postProcessing/residuals/0/residuals.dat
#gdje su dodatne oznake:
#-l oznaka za logaritamski y skalu (koristi se za residuale)
#-r 2 refreshanje grafa svake 2 sekunde
#-i 5 možeš ga zatvorit 5 sekundi nakon što se prestanu zapisivat podaci

# funkcija simulacija prima vektor brojeva koju rasporeduje po 360 stupnjeva u rx i ry
# i broj posto ce se vrtit paralelno
# kopira se cijeli folder most pod novim nazivom most1 ili most2 i td ovisno o fileNumberu tako da folder most ostaje uvijek netaknut
# zatim se izacunaju karakteristicne tocke mreze tj. tocke blokova koje se zatim upisuju u file blockMeshDict
# otvara se controlDict u koji se upisuje l_ref i A_ref za izracunavanje sila
# brisu se folderi postProcessing i dynamicCode ako postoje
# zapocinje simulacija koja racuna lift i drag force
# nakon proracuna se brise folder u kojem se izvodila simulacija
import time
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

def simulacija(X,fileNumber):

    newFileName = './proracun/proracun/most_'+str(fileNumber)
    shutil.copytree('./proracun/best/best',newFileName)
    #print(newFileName)
    #dubina uronjene noge
    h=10.
    #posto je profil simetrican racuna se samo gornja krivulja od 0 do 180 rx ide s lijeva
    n = 10      #interpolacijske tocke
    r_temp = X
    fi_temp = np.linspace(0., pi, len(r_temp))
    #print(len(fi_temp), len(r_temp))
    interp_r = interpolate.interp1d(fi_temp,r_temp,kind='cubic')
    
    
    def func_max_x(x,xt):
        delta_x=[]
        for i in range(len(x)):
            delta_x.append(abs(x[i]-xt))
        return max(delta_x)
    
    def r_phi(phi):
        phi=np.abs(phi)
        phi=phi%(2*pi)
        if(phi<pi):
            return interp_r(phi)
        else:
            return interp_r(2*pi-phi)

    def ogranicenje_superponirano (qx, qy, r, fi, Ix, Iy):
        M_x = qy * h**2 /2.
        M_y = qx * h**2 /2.
#        M_x = 10000
#        M_y = 10000
        
        x=[]
        y=[]
        x_prvi_kvadr=[]
        y_prvi_kvadr=[]
        x_drugi_kvadr=[]
        y_drugi_kvadr=[]
        x_treci_kvadr=[]
        y_treci_kvadr=[]
        x_cetvrti_kvadr=[]
        y_cetvrti_kvadr=[]
        for i in range (len(fi)):
            if (fi[i]<=pi/2):
                x_prvi_kvadr.append(r[i]*cos(fi[i]))
                y_prvi_kvadr.append(r[i]*sin(fi[i]))
                x.append(x_prvi_kvadr[-1])
                y.append(y_prvi_kvadr[-1])
            elif (fi[i]>=pi/2 and fi[i]<pi):
                x_drugi_kvadr.append(r[i]*cos(fi[i]))
                y_drugi_kvadr.append(r[i]*sin(fi[i]))
                x.append(x_drugi_kvadr[-1])
                y.append(y_drugi_kvadr[-1])
            elif (fi[i]>=pi and fi[i]<3*pi/2):
                x_treci_kvadr.append(r[i]*cos(fi[i]))
                y_treci_kvadr.append(r[i]*sin(fi[i]))
                x.append(x_treci_kvadr[-1])
                y.append(y_treci_kvadr[-1])
            else:
                x_cetvrti_kvadr.append(r[i]*cos(fi[i]))
                y_cetvrti_kvadr.append(r[i]*sin(fi[i]))
                x.append(x_cetvrti_kvadr[-1])
                y.append(y_cetvrti_kvadr[-1])
        
        sigma_x=[]
        sigma_y=[]
        for i in range (len(x_prvi_kvadr)):
            sigma_x.append(- M_x / Ix * abs(y_prvi_kvadr[i]))
            sigma_y.append(- M_y / Iy * abs(x_prvi_kvadr[i]))
        for i in range (len(x_drugi_kvadr)):
            sigma_x.append(- M_x / Ix * abs(y_drugi_kvadr[i]))
            sigma_y.append( M_y / Iy * abs(x_drugi_kvadr[i]))  
        for i in range (len(x_treci_kvadr)):
            sigma_x.append( M_x / Ix * abs(y_treci_kvadr[i]))
            sigma_y.append( M_y / Iy * abs(x_treci_kvadr[i]))
        for i in range (len(x_cetvrti_kvadr)):
            sigma_x.append( M_x / Ix * abs(y_cetvrti_kvadr[i]))
            sigma_y.append(- M_y / Iy * abs(x_cetvrti_kvadr[i]))
        
        sigma_ekv_plus=[]
        sigma_ekv_minus=[]
        for i in range(len(sigma_x)):
            sigma_ekv_plus.append(sigma_x[i]+sigma_y[i])
            sigma_ekv_minus.append(-sigma_x[i]+sigma_y[i])
        
        sigma_ekv_plus_max=max([abs(max(sigma_ekv_plus)),abs(min(sigma_ekv_plus))])
        sigma_ekv_minus_max=max([abs(max(sigma_ekv_minus)),abs(min(sigma_ekv_minus))])
        sigma_superponirano =max([sigma_ekv_minus_max,sigma_ekv_plus_max])
#        print(sigma_ekv_plus_max)
#        print(sigma_ekv_minus_max)
        
        
#        plt.plot(x,y)
#        fi_for_plot=np.rad2deg(fi)
#        plt.plot(fi_for_plot, sigma_x, label='sigma_lift')
#        plt.plot(fi_for_plot, sigma_y, label='sigma_drag')
#        plt.plot(fi_for_plot, sigma_ekv_plus, label='sigma_superpon_plus')
#        plt.plot(fi_for_plot, sigma_ekv_minus, label='sigma_superpon_minus')
#        plt.xlabel('fi')
#        plt.ylabel('Sigma')
#        plt.legend()
#        plt.savefig('sigme.png') 
        return sigma_superponirano*20

    def ogranicenje (qx, qy, x_max, y_max, Ix, Iy):
        Mx = qy * h**2 /2.
        My = qx * h**2 /2.
        sigma_x = Mx * y_max/Ix
        sigma_y = My * x_max/Iy
        return sigma_x , sigma_y
    
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
    
    fi_prvi = np.linspace(0.,pi/4, n/2+1)
    fi_drugi = np.linspace(pi/4, 3*pi/4, n+2)
    fi_treci = np.linspace(3*pi/4, pi, n/2+1)
    r_prvi = interp_r(fi_prvi[0:-1])
    x_prvi =[]
    y_prvi =[]
    for i in range (len(r_prvi)):
        x_prvi.append(r_prvi[i]*cos(fi_prvi[i]))
        y_prvi.append(r_prvi[i]*sin(fi_prvi[i]))
        #print(np.rad2deg(fi_prvi[i]))
    r_drugi = interp_r(fi_drugi[1:-1])
    x_drugi = []
    y_drugi = []
    for i in range (len(r_drugi)):
        x_drugi.append(r_drugi[i]*cos(fi_drugi[i+1]))
        y_drugi.append(r_drugi[i]*sin(fi_drugi[i+1]))
        #print(np.rad2deg(fi_drugi[i+1]))
    r_treci = interp_r(fi_treci[1:])
    x_treci =[]
    y_treci =[]
    for i in range (len(r_treci)):
        x_treci.append(r_treci[i]*cos(fi_treci[i+1]))
        y_treci.append(r_treci[i]*sin(fi_treci[i+1]))
        #print(np.rad2deg(fi_treci[i+1]))
    
    # Moment inercije Ix, Iy
    fi_fi=np.linspace(0.,2*pi,101)
    R=[]
    for i in range (len(fi_fi)):
        R.append((r_phi(fi_fi[i]))**2/2)
    Area=integrate.simps(R,fi_fi)        
    print(Area)

    def integrand_x(fi_kut):
        return 1/4.*(r_phi(fi_kut))**4*(sin(fi_kut))**2
        
    r_x=[]
    for i in range(len(fi_fi)):
        r_x.append(integrand_x(fi_fi[i]))
    Ix = integrate.simps(r_x,fi_fi)
    #print(Ix)
    ############################
#    def fr_y(r,phi):
#        return r**3*(cos(phi))**2
    #Iy,err2=integrate.nquad(fr_y,[bound_r,bound_phi])
    #print(Iy,err2)
    #provjera
    ##########################
    def integrand_y(fi_kut):
        return 1/4.*(r_phi(fi_kut))**4*(cos(fi_kut))**2
    r_y=[]
    for i in range(len(fi_fi)):
        r_y.append(integrand_y(fi_fi[i]))
    Iy = integrate.simps(r_y,fi_fi)
    #print(Iy)
    ##########################
#    def fr_y(r,phi):
#        return r*sin(phi)
    #Ty,err3=integrate.nquad(fr_y,[bound_r,bound_phi])
    #print(Ty)
    #provjera
    ##########################
    def integrand_Ty(fi_kut):
        return 1/3.*(r_phi(fi_kut))**3*sin(fi_kut)
    r_Ty=[]
    for i in range(len(fi_fi)):
        r_Ty.append(integrand_Ty(fi_fi[i]))
    Ty =1/Area*integrate.simps(r_Ty,fi_fi)
    print(Ty)
    ##########################
#    def fr_x(r,phi):
#        return r*cos(phi)
    #Tx,err4=integrate.nquad(fr_x,[bound_r,bound_phi])
    #print(Tx)
    #provjera
    ##########################
    def integrand_Tx(fi_kut):
        return 1/3.*(r_phi(fi_kut))**3*cos(fi_kut)
    r_Tx=[]
    for i in range(len(fi_fi)):
        r_Tx.append(integrand_Tx(fi_fi[i]))
    Tx = 1/Area*integrate.simps(r_Tx,fi_fi)
    print(Tx)
    #print('Tx='+str(Tx)+' Ty='+str(Ty))
    
    #domena simulacije noge mosta
    domena_x = [-45., 115.]
    domena_y = [-50., 50.]
    poddomena_x = [-15., 15.]
    poddomena_y = [-15., 15.]
    podjela_poddomena_x=30
    podjela_poddomena_y=50
    podjela_domena_x_lijevo=20
    podjela_domena_x_desno=50
    podjela_domena_y=20
    gradacija=20
    
    #trazenje karakteristinih tocki za blockMesh
    prva=[interp_r(pi/4)*cos(pi/4),interp_r(pi/4)*sin(pi/4)]
    druga=[interp_r(3*pi/4)*cos(3*pi/4),interp_r(3*pi/4)*sin(3*pi/4)]
    treca=[druga[0],-druga[1]]
    cetvrta=[prva[0],-prva[1]]
    
    #interpolacijske tocme izmedu glavnih tocaka
    interp_prva_cetvrta=[]
    for i in reversed(range(len(r_prvi))):
        interp_prva_cetvrta.append([x_prvi[i],y_prvi[i],-0.5])
    for i in range(len(r_prvi)-1):
        interp_prva_cetvrta.append([x_prvi[i+1],-y_prvi[i+1],-0.5])
    interp_prva_druga=[]
    for i in range (len(r_drugi)):
        interp_prva_druga.append([x_drugi[i],y_drugi[i],-0.5])        
    interp_druga_treca=[]
    for i in range(len(r_treci)):
        interp_druga_treca.append([x_treci[i],y_treci[i],-0.5])
    for i in reversed (range(len(interp_druga_treca)-1)):
        interp_druga_treca.append([interp_druga_treca[i][0],-interp_druga_treca[i][1],-0.5])
    interp_treca_cetvrta=[]
    for i in reversed(range(len(interp_prva_druga))):
        interp_treca_cetvrta.append([interp_prva_druga[i][0],-interp_prva_druga[i][1],-0.5])
   ####################max_x,max_y
    xevi=[]
    yloni=[]
    for i in range(len(interp_prva_cetvrta)):
        xevi.append(interp_prva_cetvrta[i][0])
        yloni.append(interp_prva_cetvrta[i][1])
    for i in range(len(interp_prva_druga)):
        xevi.append(interp_prva_druga[i][0])
        yloni.append(interp_prva_druga[i][1])
    for i in range(len(interp_druga_treca)):
        xevi.append(interp_druga_treca[i][0])
        yloni.append(interp_druga_treca[i][1])
    for i in range(len(interp_treca_cetvrta)):
        xevi.append(interp_treca_cetvrta[i][0])
        yloni.append(interp_treca_cetvrta[i][1])
    x_max = func_max_x(xevi,Tx)
    y_max = max(yloni)
#    print ('x_max= '+str(x_max)+' ,y_max='+str(y_max))
############################################################
#Pisanje blockMesh-a
############################################################
    
    points = []
    points.append(prva)                                                      #0
    points.append(druga)                                                     #1
    points.append(treca)                                                     #2
    points.append(cetvrta)                                                   #3
    points.append([poddomena_x[1],poddomena_y[1]])                           #4
    points.append([poddomena_x[0],poddomena_y[1]])                           #5
    points.append([poddomena_x[0],poddomena_y[0]])                           #6
    points.append([poddomena_x[1],poddomena_y[0]])                           #7
    points.append([domena_x[1],domena_y[1]])                                 #8
    points.append([poddomena_x[1],domena_y[1]])                              #9
    points.append([poddomena_x[0],domena_y[1]])                              #10
    points.append([domena_x[0],domena_y[1]])                                 #11
    points.append([domena_x[0],poddomena_y[1]])                              #12
    points.append([domena_x[0],poddomena_y[0]])                              #13
    points.append([domena_x[0],domena_y[0]])                                 #14
    points.append([poddomena_x[0],domena_y[0]])                              #15
    points.append([poddomena_x[1],domena_y[0]])                              #16
    points.append([domena_x[1],domena_y[0]])                                 #17
    points.append([domena_x[1],poddomena_y[0]])                              #18
    points.append([domena_x[1],poddomena_y[1]])                              #19

    spline=[]
    spline.append("    spline 0 3 (")
    spline.append("    spline 0 1 (")
    spline.append("    spline 1 2 (")
    spline.append("    spline 2 3 (")
    spline.append("    spline 20 23 (")
    spline.append("    spline 20 21 (")
    spline.append("    spline 21 22 (")
    spline.append("    spline 22 23 (")
#    print(interp_prva_cetvrta[0])
    for i in range(len(interp_prva_cetvrta)):
        spline[0]=spline[0]+"("+str(interp_prva_cetvrta[i][0])+" "+str(interp_prva_cetvrta[i][1])+" "+"-0.5) "
        spline[4]=spline[4]+"("+str(interp_prva_cetvrta[i][0])+" "+str(interp_prva_cetvrta[i][1])+" "+"0.5) "
    spline[0]=spline[0][0:-1]+")"
    spline[4]=spline[4][0:-1]+")"
    for i in range(len(interp_prva_druga)):
        spline[1]=spline[1]+"("+str(interp_prva_druga[i][0])+" "+str(interp_prva_druga[i][1])+" "+"-0.5) "
        spline[5]=spline[5]+"("+str(interp_prva_druga[i][0])+" "+str(interp_prva_druga[i][1])+" "+"0.5) "
    spline[1]=spline[1][0:-1]+")"
    spline[5]=spline[5][0:-1]+")"
    for i in range(len(interp_druga_treca)):
        spline[2]=spline[2]+"("+str(interp_druga_treca[i][0])+" "+str(interp_druga_treca[i][1])+" "+"-0.5) "
        spline[6]=spline[6]+"("+str(interp_druga_treca[i][0])+" "+str(interp_druga_treca[i][1])+" "+"0.5) "
    spline[2]=spline[2][0:-1]+")"
    spline[6]=spline[6][0:-1]+")"
    for i in range(len(interp_treca_cetvrta)):
        spline[3]=spline[3]+"("+str(interp_treca_cetvrta[i][0])+" "+str(interp_treca_cetvrta[i][1])+" "+"-0.5) "
        spline[7]=spline[7]+"("+str(interp_treca_cetvrta[i][0])+" "+str(interp_treca_cetvrta[i][1])+" "+"0.5) "
    spline[3]=spline[3][0:-1]+")"
    spline[7]=spline[7][0:-1]+")"
#    print(spline[3])
    
    target = open('./'+newFileName+'/system/blockMeshDict', 'r')
    contents = target.readlines()
    target.close()
    for j in range (len(contents)):
        # pisanje tocaka
        if contents [j] == '        pointField points(20);\n':
            for i in range (len (points)):
                contents [j+i+1] = "        points[" + str(i) + "]  = point(" + str(points[i][0]) +", "+ str(points[i][1]) + ", -0.5);\n"
        if contents [j] =='edges\n':
            for i in range (len(spline)):
                contents[j+i+2] =spline[i]+"\n"
        if contents [j] =='blocks\n':
            contents [j+2] ='    hex (3 7 4 0 23 27 24 20) ('+str(podjela_poddomena_y)+' '+str(podjela_poddomena_x)+' 1) edgeGrading ('+str(gradacija)+' 1 1)\n'
            contents [j+3] ='    hex (0 4 5 1 20 24 25 21) ('+str(podjela_poddomena_y)+' '+str(podjela_poddomena_x)+' 1) edgeGrading ('+str(gradacija)+' 1 1)\n'
            contents [j+4] ='    hex (2 1 5 6 22 21 25 26) ('+str(podjela_poddomena_x)+' '+str(podjela_poddomena_y)+' 1) edgeGrading (1 '+str(gradacija)+' 1)\n'
            contents [j+5] ='    hex (2 6 7 3 22 26 27 23 22) ('+str(podjela_poddomena_y)+' '+str(podjela_poddomena_x)+' 1) edgeGrading ('+str(gradacija)+' 1 1)\n'
            contents [j+6] ='    hex (12 5 10 11 32 25 30 31) ('+str(podjela_domena_x_lijevo)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            contents [j+7] ='    hex (5 4 9 10 25 24 29 30) ('+str(podjela_poddomena_x)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            contents [j+8] ='    hex (4 19 8 9 24 39 28 29) ('+str(podjela_domena_x_desno)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            contents [j+9] ='    hex (13 6 5 12 33 26 25 32) ('+str(podjela_domena_x_lijevo)+' '+str(podjela_poddomena_x)+' 1) simpleGrading (1 1 1)\n'
            contents [j+10] ='    hex (7 18 19 4 27 38 39 24) ('+str(podjela_domena_x_desno)+' '+str(podjela_poddomena_x)+' 1) simpleGrading (1 1 1)\n'
            contents [j+11] ='    hex (14 15 6 13 34 35 26 33) ('+str(podjela_domena_x_lijevo)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            contents [j+12] ='    hex (15 16 7 6 35 36 27 26) ('+str(podjela_poddomena_x)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            contents [j+13] ='    hex (16 17 18 7 36 37 38 27) ('+str(podjela_domena_x_desno)+' '+str(podjela_domena_y)+' 1) simpleGrading (1 1 1)\n'
            #print(contents [j+13])
        
    out = open('./'+newFileName+'/system/blockMeshDict', 'w')
    for i in range(len(contents)):
        out.writelines(str(contents[i]))
    out.close()
    
    l_ref=abs(X[0])+abs(X[-1])
    ry=[]
    for i in range (len(r_temp)):
        ry.append(r_temp[i]*sin(fi_temp[i]))
    A_ref=max(ry)*2.
    
    target = open('./'+newFileName+'/system/controlDict', 'r')
    contents_controlDict = target.readlines()
    target.close()
    contents_controlDict[73] = '        lRef        '+str(l_ref)+';\n'
    contents_controlDict[74] = '        ARef        '+str(A_ref)+';\n'
    
    out_controlDict = open('./'+newFileName+'/system/controlDict', 'w')
    for i in range(len(contents_controlDict)):
        out_controlDict.writelines(str(contents_controlDict [i]))
    out_controlDict.close()

    case_dir = newFileName
    solver = 'icoFoam'
    
    """
    Start OpenFOAM simulation.
    """
    
    
    print(' - Starting calculation... %s' % case_dir)
    
    #brisanje foldera postProcessing i dynamicCode ako postoji
    list_file = os.listdir('./'+newFileName)
    for i in range (len(list_file)):
        if list_file[i] == 'postProcessing':
            subprocess.call(['rm', '-rf', './'+newFileName+'/postProcessing'])
        elif list_file[i] == 'dynamicCode':
            subprocess.call(['rm', '-rf', './'+newFileName+'/dynamicCode'])
    
    cmd = 'bash start_case.sh %s %s' % (case_dir, solver)
    proc = subprocess.Popen(cmd.split(), cwd='./', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print (out, err)
    print(' - Calculation finished! %s' % case_dir)
    exit_code = proc.wait()
    print(' - Return code:', exit_code)
    
    if exit_code != 0:
        #shutil.rmtree(newFileName)
        return 1e15, 1e15, 1e10*Area, 1e15, 1e15, 1e10
    
    forces_file = './'+newFileName+'/postProcessing/forcesIncompressible/0/forces.dat'
    
    if not os.path.isfile(forces_file):
        print ("Forces file not found at "+forces_file)
        print ("Be sure that the case has been run and you have the right directory!")
        print ("Exiting.")
        #shutil.rmtree(newFileName)
        sys.exit()
    
    time = []
    drag = []
    lift = []
    moment = []
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
    #plt.plot(time[5:],lift[5:])
    
    #   maksimalne sile po metru visine mosta    
    #drag_max = (max([max(drag[500:-1]),abs(min(drag[500:-1]))]))+1000
    #lift_max = (max([max(lift[500:-1]),abs(min(lift[500:-1]))]))+1000
    
    drag_max = np.max(np.abs(drag[500:])) + 1000.
    lift_max = np.max(np.abs(lift[500:])) + 1000.
    
    n_points= 361
    fi=np.linspace(0.,2*pi,n_points)
    r=[]
    for i in range (len(fi)):
        r.append(r_phi(fi[i]))
#    plt.plot(x_prvi_kvadr,y_prvi_kvadr)
#    plt.plot(x_drugi_kvadr,y_drugi_kvadr)           
    
    Sigma_max = ogranicenje_superponirano (drag_max, lift_max, r, fi, Ix, Iy)
    #Sigma_x, Sigma_y = ogranicenje (drag_max, lift_max, x_max, y_max, Ix, Iy)          

    #shutil.rmtree(newFileName)
    #print(Sigma_x,Sigma_y)
    #print (Sigma_x/10**5, Sigma_y/10**5, Area)
    return Sigma_max, Area, drag_max, lift_max, Ix, Iy

faktor_sigurnosti =1.
sigma_dop = 2*10**5/faktor_sigurnosti

if __name__ == "__main__":
    def kvadr_pop_presjek (b,h):
        fi_diag=math.atan2(h/2,b/2)
        #print(np.rad2deg(fi_diag))
        fi=np.linspace(0,pi,360)
        r=[]
        for i in range (len(fi)):
            if (fi[i]<=fi_diag):
                x=b/2
                r.append(x/np.cos(fi[i]))
            elif (fi[i]>fi_diag and fi[i]<=pi-fi_diag):
                y=h/2
                r.append(y/np.sin(fi[i]))
            else:
                x=b/2
                r.append(abs(x/np.cos(fi[i])))
        return r
    #X=[ 2., 2.2, 1.8, 1.8, 2.,2.,2.,2.]
    #X=[3.69803155,3.78795551,2.09015827,4.44300334,1.08248567,3.67183216,0.59487738,3.24516277]
    #X=[4.38561015,2.11143417,5.3934392,4.44132813,2.04115363,5.51451399,3.67158115,1.35242348]
    #X=[2.50690218,4.92808389,2.86215534,1.9794298,5.63875802,4.49786406,4.33377154,1.7460372]
    #X=[3.69567485,1.11910148,4.08031717,2.60037542,1.05344528,2.62814699,5.04284431,1.58244596]
    #X=[2,2,2,2,2,2,2,2]
    X=[3.2207603893825074, 0.90801560199401887, 3.157813375761541, 2.7724328583924893, 2.0377827745169523, 1.8084522673720427, 1.9963162071726572, 6.0]
    start = time.time()
    #X= kvadr_pop_presjek(2,4)
    sigma_max, A, drag_max, lift_max, Ix, Iy =simulacija(X,'gg')
    end = time.time() 
    print(sigma_max, A, drag_max, lift_max, Ix, Iy)
    print (end-start)

#sigma_x, sigma_y, cilj = simulacija(X,1)
#print ((sigma_dop -sigma_x)*10**-6, (sigma_dop - sigma_y)*10**-6, cilj)


