# -- coding: utf-8 --
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import subprocess as sp
import os
from StringIO import StringIO
import array
import random
from scipy import stats

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import pylab as pl
import matplotlib.pyplot as plt
import scipy.interpolate
from math import pi, atan, acos, cos, sin, tan, exp, sqrt, floor
from siteOpt import Airfoil, readDAT, AoAsteps, BlendAirfoils, Turbine, individualNr, readTUR

# Python v 2.7

# Initiate fixed variables

R = 10.046/2
hubRadius = 0.72 # radius where the blade begins and hub ends (m)

TSR = 8 # tip speed ratio (a.k.a. the symbol lambda)
RPM = 71.63 # set only if constant RPM is to be used! otherwise set to None
B = 3 # number of blades
visc = 1.5e-5 # kinematic viscosity of air 
rho = 1.225 # density of air 
tol = 1.e-4 # convergence tolerance for a and adash
cutIn = 5
cutOut = 25
ratedPower = 19800 # max generator power (W)

chord = [0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572]
twist = [44, 44, 36, 24, 20, 16, 13, 10, 9, 4, 2, 1, 0] 
pitch = 3

# Wind distribution
windSpeeds = np.array([0.500, 1.500, 2.500, 3.500, 4.500, 5.500, 6.500, 7.500, 8.500, 9.500, 10.500, 11.500, 12.500, 13.500, 14.500, 15.500, 16.500, 17.500, 18.500, 19.500, 20.500, 21.500, 22.500])
windObs = np.array([0.009, 0.038, 0.115, 0.105, 0.111, 0.117, 0.090, 0.083, 0.080, 0.056, 0.048, 0.040, 0.026, 0.025, 0.016, 0.013, 0.008, 0.006, 0.004, 0.002, 0.001, 6.570e-4, 1.349e-4])
windFrq = windObs / np.sum(windObs) # Normalized so that area under = 1

cordMax = max(chord)
cordMin = min(chord)
Remax = sqrt(cutOut**2 + (TSR*cutOut)**2)*cordMax/visc
Remin = sqrt(cutIn**2 + (TSR*cutIn)**2)*cordMin/visc

print "remax", Remax/100000
print "remin", Remin/100000

#Re = np.linspace(0.1*Remin/100000, Remax/100000, 3) # At wich Re*10^5 to evaluate. The more the better.
Re = np.array([20])




# Set up GA with package DEAP for a single objective optimization

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

# Function evaluating rotor anunal yield

def fitnessFunc(individual):

    individualNumber = str(individualNr.next()) # Counter used in naming files with individual number

    # Split up the genome in nice chunks
    twistAndChord = individual[66:]
    twistData = twistAndChord[:4]
    chordData = twistAndChord[4:]
    foilData = individual[:66]
    foilData = np.split(foilData, 3) 
    
    rootFoil = foilData[0]
    rootFoil = np.split(rootFoil, 2) 
    rootFoilx = rootFoil[0]
    rootFoily = rootFoil[1]

    midFoil = foilData[1]
    midFoil = np.split(midFoil, 2) 
    midFoilx = midFoil[0]
    midFoily = midFoil[1]

    topFoil = foilData[2]
    topFoil = np.split(topFoil, 2) 
    topFoilx = topFoil[0]
    topFoily = topFoil[1]


    try:
        AFroot = Airfoil(rootFoilx, rootFoily, Re=Re)
        AFmid = AFroot
        AFtop = AFroot

        top = AFroot.y[3:(floor(len(AFroot.y)/2)-3)]
        bott = AFroot.y[(4+floor(len(AFroot.y)/2)):-3]
        bott = bott[::-1]

        # If a y-coordinate that should be above is below
        if False in (top > bott): 
            return [(-99999)]
        """
        fig = pl.figure(figsize=(9, 3.4))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'x/c')
        ax.set_ylabel(r'y/c')
        pl.plot(AFroot.originalx, AFroot.originaly, 'ko')
        pl.plot(AFroot.x, AFroot.y, 'r')
        plt.plot(1, 0, 'ro')
        x1,x2,y1,y2 = pl.axis()
        pl.axis((x1,x2,y1,y2))
        pl.grid()
        pl.show()
        fig.savefig('connectTheDots3.eps', dpi=fig.dpi)
        """

        AFmid = AFroot
        AFtop = AFroot
        #AFmid = Airfoil(midFoilx, midFoily, Re=Re)
        #AFtop = Airfoil(topFoilx, topFoily, Re=Re)
    except:
        return [(-99999)] # if xfoil fails

    """
    theTurbine = Turbine([AFroot, AFmid, AFtop], R=R, hubRadius=hubRadius, TSR=TSR, 
                          ratedPower=ratedPower, B=B, visc=visc, rho=rho, tol=tol, 
                          windSpeeds=windSpeeds, windFrq=windFrq, cutIn=cutIn, 
                          cutOut=cutOut, skipBlending=1, RPM=RPM, metaData=individual, indNr=individualNumber,
                          chord=chord, twist=twist, pitch=pitch)

    """
    # Use twist and chorddist from GA
    theTurbine = Turbine([AFroot, AFmid, AFtop], R=R, hubRadius=hubRadius, TSR=TSR, 
                          ratedPower=ratedPower, B=B, visc=visc, rho=rho, tol=tol, 
                          windSpeeds=windSpeeds, windFrq=windFrq, cutIn=cutIn, 
                          cutOut=cutOut, skipBlending=1, RPM=RPM, metaData=individual, indNr=individualNumber,
                          chord=chord, twist=twist, pitch=pitch,
                          twistData=twistData.tolist(), chordData=chordData.tolist())
    
    
    try:
        #print "Turbine avg power: " + str(theTurbine.avgPower()) + " W, Ind no: " + individualNumber + ". " + str(theTurbine.power)
        print "Turbine avg power: " + str(theTurbine.avgPower()) + " W, Ind no: " + individualNumber 

        return [(theTurbine.avgPower())]
    except:
        return [(-99999)]


def addPercent(value, percent):
    return value*(percent/100 + 1)

def removePercent(value, percent):
    return value*(1 - percent/100)

diff = 10 # percent

# Generate attributes around a S809 - Root

# Top side
toolbox.register("Rx2", random.uniform, removePercent(0.801, diff), addPercent(0.801, diff))
toolbox.register("Ry2", random.uniform, removePercent(0.038, diff), addPercent(0.038, diff))
toolbox.register("Rx3", random.uniform, removePercent(0.647, diff), addPercent(0.647, diff))
toolbox.register("Ry3", random.uniform, removePercent(0.068, diff), addPercent(0.068, diff))
toolbox.register("Rx4", random.uniform, removePercent(0.413, diff), addPercent(0.413, diff))
toolbox.register("Ry4", random.uniform, removePercent(0.101, diff), addPercent(0.101, diff))
toolbox.register("Rx5", random.uniform, removePercent(0.236, diff), addPercent(0.236, diff))
toolbox.register("Ry5", random.uniform, removePercent(0.089, diff), addPercent(0.089, diff))
toolbox.register("Rx6", random.uniform, removePercent(0.1, diff), addPercent(0.1, diff))
toolbox.register("Ry6", random.uniform, removePercent(0.059, diff), addPercent(0.059, diff))

# Front
toolbox.register("Rx7", random.uniform, removePercent(0, 0), addPercent(0, 0))
toolbox.register("Ry7", random.uniform, removePercent(0, 0), addPercent(0, 0))

# Bottom side
toolbox.register("Rx8", random.uniform, removePercent(0.1, diff), addPercent(0.1, diff))
toolbox.register("Ry8", random.uniform, removePercent(-0.056, diff), addPercent(-0.05, diff))
toolbox.register("Rx9", random.uniform, removePercent(0.236, diff), addPercent(0.2, diff))
toolbox.register("Ry9", random.uniform, removePercent(-0.094, diff), addPercent(-0.09, diff))
toolbox.register("Rx10", random.uniform, removePercent(0.413, diff), addPercent(0.3, diff))
toolbox.register("Ry10", random.uniform, removePercent(-0.107, diff), addPercent(-0.11, diff))
toolbox.register("Rx11", random.uniform, removePercent(0.647, diff), addPercent(0.5, diff)) 
toolbox.register("Ry11", random.uniform, removePercent(-0.056, diff), addPercent(-0.094, diff))
toolbox.register("Rx12", random.uniform, removePercent(0.801, diff), addPercent(0.8, diff))
toolbox.register("Ry12", random.uniform, removePercent(-0.020, diff), addPercent(-0.02, diff))

# Generate attributes around a naca4418 - Mid

# Top side
toolbox.register("Mx2", random.uniform, 0.7, 0.9)
toolbox.register("My2", random.uniform, 0.06, 0.07)
toolbox.register("Mx3", random.uniform, 0.4, 0.6)
toolbox.register("My3", random.uniform, 0.1, 0.15)
toolbox.register("Mx4", random.uniform, 0.2, 0.3)
toolbox.register("My4", random.uniform, 0.10, 0.15)
toolbox.register("Mx5", random.uniform, 0.09, 0.11)
toolbox.register("My5", random.uniform, 0.085, 0.095)
toolbox.register("Mx6", random.uniform, 0.023, 0.027) 
toolbox.register("My6", random.uniform, 0.048, 0.052)

# Front
toolbox.register("Mx7", random.uniform, 0.0001, 0.0002) 
toolbox.register("My7", random.uniform, 0.0001, 0.0002)

# Bottom side
toolbox.register("Mx8", random.uniform, 0.012, 0.0128) 
toolbox.register("My8", random.uniform, -0.0001, -0.020)
toolbox.register("Mx9", random.uniform, 0.0730, 0.0770) 
toolbox.register("My9", random.uniform, -0.02, -0.05)
toolbox.register("Mx10", random.uniform, 0.19, 0.21) 
toolbox.register("My10", random.uniform, -0.054, -0.056)
toolbox.register("Mx11", random.uniform, 0.39, 0.41) 
toolbox.register("My11", random.uniform, -0.056, -0.048)
toolbox.register("Mx12", random.uniform, 0.69, 0.71) 
toolbox.register("My12", random.uniform, -0.048, -0.025)

# Generate attributes around a naca4418 - Top

# Top side
toolbox.register("Tx2", random.uniform, 0.7, 0.9)
toolbox.register("Ty2", random.uniform, 0.06, 0.07)
toolbox.register("Tx3", random.uniform, 0.4, 0.6)
toolbox.register("Ty3", random.uniform, 0.1, 0.15)
toolbox.register("Tx4", random.uniform, 0.2, 0.3)
toolbox.register("Ty4", random.uniform, 0.10, 0.15)
toolbox.register("Tx5", random.uniform, 0.09, 0.11)
toolbox.register("Ty5", random.uniform, 0.085, 0.095)
toolbox.register("Tx6", random.uniform, 0.023, 0.027) 
toolbox.register("Ty6", random.uniform, 0.048, 0.052)

# Front
toolbox.register("Tx7", random.uniform, 0.0001, 0.0002) 
toolbox.register("Ty7", random.uniform, 0.0001, 0.0002)

# Bottom side
toolbox.register("Tx8", random.uniform, 0.012, 0.0128) 
toolbox.register("Ty8", random.uniform, -0.0001, -0.020)
toolbox.register("Tx9", random.uniform, 0.0730, 0.0770) 
toolbox.register("Ty9", random.uniform, -0.02, -0.05)
toolbox.register("Tx10", random.uniform, 0.19, 0.21) 
toolbox.register("Ty10", random.uniform, -0.054, -0.056)
toolbox.register("Tx11", random.uniform, 0.39, 0.41) 
toolbox.register("Ty11", random.uniform, -0.056, -0.048)
toolbox.register("Tx12", random.uniform, 0.69, 0.71) 
toolbox.register("Ty12", random.uniform, -0.048, -0.025)

# Twist 
toolbox.register("Tmax", random.uniform, 45, 50)
toolbox.register("Tx", random.uniform, 0.7, 0.9)
toolbox.register("Ty", random.uniform, 10, 2)
toolbox.register("Tmin", random.uniform, 1, -1)

# Chord 
toolbox.register("Cmax", random.uniform, 0.40, 0.42)
toolbox.register("Cx", random.uniform, 0.4, 0.6)
toolbox.register("Cy", random.uniform, 0.40, 0.42)
toolbox.register("Cmin", random.uniform, 0.40, 0.42)



toolbox.register("individual", tools.initCycle, creator.Individual,
                 (toolbox.Rx2, toolbox.Rx3, toolbox.Rx4, toolbox.Rx5, toolbox.Rx6, toolbox.Rx7, 
                  toolbox.Rx8, toolbox.Rx9, toolbox.Rx10, toolbox.Rx11, toolbox.Rx12, 
                  toolbox.Ry2, toolbox.Ry3, toolbox.Ry4, toolbox.Ry5, toolbox.Ry6, toolbox.Ry7, 
                  toolbox.Ry8, toolbox.Ry9, toolbox.Ry10, toolbox.Ry11, toolbox.Ry12, 

                  toolbox.Mx2, toolbox.Mx3, toolbox.Mx4, toolbox.Mx5, toolbox.Mx6, toolbox.Mx7, 
                  toolbox.Mx8, toolbox.Mx9, toolbox.Mx10, toolbox.Mx11, toolbox.Mx12, 
                  toolbox.My2, toolbox.My3, toolbox.My4, toolbox.My5, toolbox.My6, toolbox.My7, 
                  toolbox.My8, toolbox.My9, toolbox.My10, toolbox.My11, toolbox.My12, 

                  toolbox.Tx2, toolbox.Tx3, toolbox.Tx4, toolbox.Tx5, toolbox.Tx6, toolbox.Tx7, 
                  toolbox.Tx8, toolbox.Tx9, toolbox.Tx10, toolbox.Tx11, toolbox.Tx12, 
                  toolbox.Ty2, toolbox.Ty3, toolbox.Ty4, toolbox.Ty5, toolbox.Ty6, toolbox.Ty7, 
                  toolbox.Ty8, toolbox.Ty9, toolbox.Ty10, toolbox.Ty11, toolbox.Ty12, 

                  toolbox.Tmax, toolbox.Tx, toolbox.Ty, toolbox.Tmin,

                  toolbox.Cmax, toolbox.Cx, toolbox.Cy, toolbox.Cmin), n=1)

toolbox.register("population", tools.initRepeat, list, toolbox.individual)

 
individualNr = individualNr()

toolbox.register("evaluate", fitnessFunc)
toolbox.register("mate", tools.cxTwoPoints)
toolbox.register("mutate", tools.mutGaussian, mu=0.0, sigma=0.02, indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)

def main():
    random.seed(63)
    
    pop = toolbox.population(n=600)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", tools.mean)
    stats.register("std", tools.std)
    stats.register("min", min)
    stats.register("max", max)
    
    algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.1, ngen=10000, stats=stats,
                        halloffame=hof, verbose=True)
    
    return pop, stats, hof

if __name__ == "__main__":
    main()




