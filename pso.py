# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:45:21 2016

@author: stefan
"""

import numpy as np
import matplotlib.pyplot as plt
import inspyred
from random import Random
from time import time
import most
from multiprocessing.pool import ThreadPool
import shutil
import os
from scipy import interpolate

number_of_threads = 4

plot_F = [] # Fitness kroz iteracije
BEST = None # Best overall

sve=[]
sve_best=[]
#F_D=[]
#F_L=[]
#sigma_x_max=[]
#sigma_y_max=[]
#I_x=[]
#I_y=[]
#A=[]

r_min=0.1
r_max=8.

donje_granice = [r_min]*8
gornje_granice = [r_max]*8

#donje_granice[0]=6
#gornje_granice[0]=6
donje_granice[-1]=6
gornje_granice[-1]=6



shutil.rmtree('./proracun/proracun')
os.mkdir('./proracun/proracun')
shutil.rmtree('./proracun/best/Best_iter')
os.mkdir('./proracun/best/Best_iter')

Random_gen_min = np.ones(len(donje_granice)) * 3.0
Random_gen_min[0]=1.
Random_gen_min[-1]=6.
Random_gen_max = np.ones(len(donje_granice)) * 4.0
Random_gen_max[0]= 4.
Random_gen_max[-1]=6.


# Random generiranje

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return most.simulacija(*a_b)

def generate(random, args):
    ndim = args.get('num_inputs') # broj opt. varijabli
    particel = np.random.uniform(Random_gen_min, Random_gen_max, ndim)
    #particle_text = open('particle.txt', 'ab')
    #particle_text.write(particel)
    #print (particel)
    return particel


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
    plt.xlim(-6, 6)
    plt.savefig(filename)
    plt.close()
    
# Funkcija cilja
def evaluate(candidates, args):
    start = time()
    pool = ThreadPool(number_of_threads)
    F = []
    global sve
    sve=[]
    
    number = []
    for i in range (len(candidates)):
        number.append(i)
    SVE_sigme_povrsine = pool.map(func_star, zip( candidates, number))
    pool.close()
    pool.join()
    
    for i in range (len(candidates)):
        p = 0.0
        A0 = 30.0
        if SVE_sigme_povrsine[i][0] > most.sigma_dop:
            p = SVE_sigme_povrsine[i][0]/most.sigma_dop * 1.0
        
        R = candidates[i]
        """
        # Uvjet glatkoće
        for ir in range(1, len(R)-1):
            p += np.max([0.5 * (R[ir-1] + R[ir+1]) - R[ir], 0.1]) - 0.1
        """      
        
        fitness = SVE_sigme_povrsine[i][1]/A0 + p
        sve.append(SVE_sigme_povrsine[i])
        F.append(fitness)
        
        profilepng = './proracun/proracun/most_%d/profile.png' % i
        #plot_profile(R, profilepng)
        
        
    end = time()
    print('time='+str(end-start))
    return F

# Funkcija za monitoring optimizacije (poziva se nakon svake iteracje)    
def observe(population, num_generations, num_evaluations, args):
    
    best = max(population) # najbolji u iteraciji
    
    global BEST, plot_F # da sigurno možemo koristiti varijable definirane izvan funkcije
    global sve, sve_best
    
    # best je najbolji u iteraciji, ne mora biti najbolji sveukupno
    if num_generations == 0 or best.fitness <= BEST.fitness: # nulta iteracija
        BEST = best
        for i in range (len(population)):
            if population[i].fitness==BEST.fitness :
                print('BEST >>>>>>>>>>> ', i)
                #print (population[i])
                #shutil.rmtree('./proracun/best/best')
                #os.makedirs('./proracun/best/best')
                best_case='./proracun/proracun/most_'+str(i)
                shutil.rmtree('./proracun/best/best/0')
                shutil.copytree(best_case+'/300','./proracun/best/best/0')
                #shutil.copytree(best_case+'/constant','./proracun/best/best/constant')
                #shutil.copytree(best_case+'/system','./proracun/best/best/system')
                shutil.copytree(best_case,'./proracun/best/Best_iter/best_'+str(num_generations))
                sve_best.append(sve[i])
                
#                shutil.copyfile('./proracun/best/Best_iter/best_%d/profile.png' % num_generations, './proracun/profile_%03d.png' % num_generations)
                
                break
    else:
        sve_best.append(sve_best[-1])
    plot_F.append(BEST.fitness)
    
    #print(sve_best)
    
    list_file = os.listdir('./proracun/proracun')
    for i in range (len(list_file)):
        shutil.rmtree('./proracun/proracun/'+list_file[i])
    #print(population[0].candidate)
    
    #upisivanje najboljih F_D,F_L,A,...
    sve_ostalo= open('./proracun/ostalo.txt', 'r')
    contents = sve_ostalo.readlines()
    sve_ostalo.close()
    contents.append(str(num_generations)+'\t'+ str(sve_best[-1][0])+'\t'+ str(sve_best[-1][1])+'\t'+ str(sve_best[-1][2])+'\t'+ str(sve_best[-1][3])+'\t'+ str(sve_best[-1][4])+'\t'+ str(sve_best[-1][5])+'\n')
    out = open('./proracun/ostalo.txt', 'w')
    for i in range(len(contents)):
        out.writelines(str(contents[i]))
    out.close()
    
    #upisivanje kandidata
    for i in range (len(population)):
        path_kandidata = './proracun/candidate_'+str(i)+'.txt'
        kandidati = open(path_kandidata, 'r')
        contents = kandidati.readlines()
        kandidati.close()
        contents.append('\n'+str(num_generations)+'\t')
        contents.append(str(population[i].candidate)+'\t'+str(population[i].fitness)+'\t'+str(sve[i][0])+'\t'+str(sve[i][1]))
        out = open(path_kandidata, 'w')
        for i in range(len(contents)):
            out.writelines(str(contents[i]))
        out.close()
        
    # crtaj ukupno najbolju putanju
#    OuT_FunkcijaCilja.f(BEST.candidate) , True, 'slike/iter_%05d.png' % num_generations)
    rezultati = open('./proracun/rezultati.txt', 'r')
    contents = rezultati.readlines()
    rezultati.close()
    contents.append(str(num_generations)+'\t'+ str(best.fitness)+'\t'+ str(BEST.fitness)+'\t'+ str(BEST.candidate)+'\n')#+ str(BEST.candidate[1])+\
                    #'\t'+ str(BEST.candidate[2])+'\t'+ str(BEST.candidate[3])+'\t'+ str(BEST.candidate[4])+'\t'+ str(BEST.candidate[5])+'\t'+\
                    #str(BEST.candidate[6])+'\t'+ str(BEST.candidate[7])+'\n')
    out = open('./proracun/rezultati.txt', 'w')
    for i in range(len(contents)):
        out.writelines(str(contents[i]))
    out.close()
    
    #print ('iter: %5d     iter_best: %10.6f     all_best: %10.6f     best_l: %10.6f     best_r1: %10.6f     best_k: %10.6f' % (num_generations, best.fitness, BEST.fitness, BEST.candidate[0], BEST.candidate[1], BEST.candidate[2]))
    
    
rand = Random()
rand.seed(time())
    
pso = inspyred.swarm.PSO(rand)
pso.terminator = inspyred.ec.terminators.evaluation_termination
pso.observer = observe

final_swarm = pso.evolve(generator = generate,
                         evaluator = evaluate,
                         maximize = False, # minimizacija
                         bounder = inspyred.ec.Bounder(donje_granice, gornje_granice), # granice varijabli (mogu se zadati liste ako su različite)
                         num_inputs = 8, # broj dimanzija (optimizacijskih varijabli)
                         pop_size = 8, # swarm size
                         max_evaluations = 6400, # maxsimalni broje evaluacije (max_iter = max_evals / swarm_size)
                         inertia = 0.7, # faktor inercije
                         cognitive_rate = 1.0, # kognitivni faktor
                         social_rate = 1.0 # socijalni faktor
                         )

final_swarm.sort(reverse=True)

#best = final_swarm[0]
#best_A = best.candidate
#best_f = best.fitness

# prikazi graf najbolje putanje
#f.f(best_x, True)

# Graf konvergencije
"""
plt.figure()
plt.plot(plot_F)
plt.xlabel('Iterations')
plt.ylabel('Fitness')
plt.show()
plt.savefig('PSO.png')
"""
