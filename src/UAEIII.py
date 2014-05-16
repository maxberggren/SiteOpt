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
from math import pi, atan, acos, cos, sin, tan, exp, sqrt
from siteOpt import Airfoil, readDAT, AoAsteps, BlendAirfoils, Turbine, individualNr, readDAT

# Python v 2.7

xs,ys = readDAT("S809.dat")

#pl.plot(xs, ys)
#pl.show()


R = 10.046/2
hubRadius = 0.72 # radius where the blade begins and hub ends (m)

TSR = 8 # tip speed ratio (a.k.a. the symbol lambda)
RPM = 71.63 # set only if constant RPM is to be used! otherwise set to None because this overrides tip speed ratio
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

Re = np.linspace(0.3*Remin/100000, Remax/100000, 12) # At wich Re*10^5 to evaluate. The more the better.
#Re = np.array([20])


AFroot = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat")
"""
AFroot = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat", 
	     alfa=[-3.09, -2.06, -1.04, -0.01, 1.02, 2.05, 3.08, 4.10, 5.13, 6.16, 7.17, 8.20, 9.22, 10.21, 11.21, 12.22, 13.24, 14.24, 15.24, 16.24, 17.23], 
	     Cldata=[-0.23, -0.10, 0.02, 0.14, 0.26, 0.39, 0.51, 0.63, 0.74, 0.82, 0.89, 0.97, 1.01, 1.00, 0.98, 1.01, 1.05, 1.08, 1.10, 1.06, 1.00],
	     Cddata=[0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.05, 0.08, 0.08, 0.09, 0.11, 0.16, 0.19])

AFroot = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat")


Cldata=[-0.066, -0.217, -0.456, -0.656, -0.767, -0.782, -0.736, -0.780, -0.793, -0.784, -0.667, -0.588, -0.578, -0.604, -0.531, -0.364, -0.211, -0.033, 0.144, 0.355, 0.783, 0.934, 0.967, 0.956, 0.951, 1.009, 1.055, 1.058, 1.036, 0.983, 0.867, 0.779, 0.751, 0.780, 0.779, 0.738, 0.648, 0.506, 0.322, 0.092, -0.127]
Cddata=[1.228, 1.206, 1.114, 0.947, 0.743, 0.528, 0.366, 0.214, 0.171, 0.076, 0.048, 0.029, 0.024, 0.019, 0.019, 0.014, 0.014, 0.005, 0.010, 0.010, 0.010, 0.019, 0.029, 0.033, 0.038, 0.052, 0.071, 0.081, 0.100, 0.124, 0.162, 0.252, 0.352, 0.513, 0.656, 0.803, 0.946, 1.079, 1.164, 1.212, 1.193]
alfa=[-88.377, -81.246, -68.344, -55.392, -43.484, -32.884, -25.337, -18.292, -16.602, -14.652, -11.774, -8.895, -8.564, -6.727, -6.108, -4.426, -3.300, -1.820, 0.214, 1.840, 5.563, 7.682, 9.116, 10.207, 11.298, 12.037, 14.050, 15.154, 16.485, 17.300, 17.920, 18.587, 24.528, 32.356, 39.113, 47.066, 55.056, 65.102, 75.287, 87.338, 98.340]

print len(Cldata), len(Cddata), len(alfa)
AFroot = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat", 
	     alfa=alfa,
	     Cldata=Cldata,
	     Cddata=Cddata,
	     highAngleInterp=0)
"""

theTurbine = Turbine([AFroot, AFroot, AFroot], R=R, hubRadius=hubRadius, TSR=TSR, 
                      ratedPower=ratedPower, B=B, visc=visc, rho=rho, tol=tol, 
                      windSpeeds=windSpeeds, windFrq=windFrq, cutIn=cutIn, 
                      cutOut=cutOut, skipBlending=1, RPM=RPM,
                      chord=chord, twist=twist, pitch=pitch)


print "Turbine avg power: " + str(theTurbine.avgPower()) 


print theTurbine.power
print theTurbine.torque
print theTurbine.RPMs
print theTurbine.Cps


fig = pl.figure(figsize=(7, 5))
ax = fig.add_subplot(111)
ax.set_xlabel(r'Vindhastighet (m/s)')
ax.set_ylabel(r'Effekt (W)')

p1, = pl.plot(np.linspace(cutIn, cutOut, 10), map(theTurbine.powerInterpolate, np.linspace(cutIn, cutOut, 10)))
p2, = pl.plot([5,6,7,8,9,10,11,12,13], [0,1500,3300,5050,7500,9500,10500,11000,11500], 'o')
p3, = pl.plot([4.975, 5.976, 6.978, 7.976, 8.986, 9.979, 10.982, 11.988, 12.976, 13.983, 14.989], [1939.900, 2548.558, 5300.329, 7065.321, 9123.519, 11986.682, 16073.783, 19063.977, 20079.064, 19843.881, 19084.537], '-')




pl.legend([p1, p2, p3], [u'Studiens modell (SiteOpt)', u'NREL experimentell data', u'BEM-data Chepyala'], loc=2)

x1,x2,y1,y2 = pl.axis()
pl.axis((x1,14,0,20000))

pl.grid()


pl.show()
fig.savefig('UAEIII.eps', dpi=fig.dpi)


fig4 = pl.figure(figsize=(8, 5))
ax5 = fig4.add_subplot(111)
p1, = pl.plot(AFroot.alfa, AFroot.Cddata[0], linestyle='-')
#p2, = pl.plot(AFS809.alfa, AFS809.Cddata[0], linestyle="--")
p3, = pl.plot([-3.09, -2.06, -1.04, -0.01, 1.02, 2.05, 3.08, 4.10, 5.13, 6.16, 7.17, 8.20, 9.22, 10.21, 11.21, 12.22, 13.24, 14.24, 15.24, 16.24, 17.23], [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.03, 0.05, 0.08, 0.08, 0.09, 0.11, 0.16, 0.19], linestyle="--", marker='o')

x1,x2,y1,y2 = pl.axis()
pl.axis((-3,150,0,2))
pl.legend([p1, p3], ["asdasdadsdas", "Experimentell data"], loc=2)
ax5.set_xlabel(r'$\alpha$')
ax5.set_ylabel(u'$C_d$')
pl.grid()
pl.show()

fig5 = pl.figure(figsize=(8, 5))
ax6 = fig5.add_subplot(111)
p1, = pl.plot(AFroot.alfa, AFroot.Cldata[0], linestyle='-')
#p2, = pl.plot(AFS809.alfa, AFS809.Cddata[0], linestyle="--")
p3, = pl.plot([-3.09, -2.06, -1.04, -0.01, 1.02, 2.05, 3.08, 4.10, 5.13, 6.16, 7.17, 8.20, 9.22, 10.21, 11.21, 12.22, 13.24, 14.24, 15.24, 16.24, 17.23], [-0.23, -0.10, 0.02, 0.14, 0.26, 0.39, 0.51, 0.63, 0.74, 0.82, 0.89, 0.97, 1.01, 1.00, 0.98, 1.01, 1.05, 1.08, 1.10, 1.06, 1.00], linestyle="--", marker='o')

x1,x2,y1,y2 = pl.axis()
pl.axis((-3,150,0,1.3))
pl.legend([p1, p3], ["asdasdadsdas", "Experimentell data"], loc=2)
ax6.set_xlabel(r'$\alpha$')
ax6.set_ylabel(u'$C_l$')
pl.grid()
pl.show()