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



#pl.plot(xs, ys)
#pl.show()


R = 3.3/2
hubRadius = 0.144 # radius where the blade begins and hub ends (m)

TSR = 5.3 # tip speed ratio (a.k.a. the symbol lambda)
RPM = 200 # set only if constant RPM is to be used! otherwise set to None because this overrides tip speed ratio
B = 3 # number of blades
visc = 1.5e-5 # kinematic viscosity of air 
rho = 1.225 # density of air 
tol = 1.e-4 # convergence tolerance for a and adash
cutIn = 3.6
cutOut = 25
ratedPower = 19800 # max generator power (W)

chord = np.linspace(0.3,0.1,13)
twist = np.linspace(18,1,13)
pitch = 0

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
#Re = np.array([5])

xs,ys = readDAT("S835.dat")
AFroot = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat")
#AFroot.Cddata[0, 60:] = 2
#print AFroot.Cddata[0]
AFmid = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat")
xs,ys = readDAT("S833.dat")
AFtop = Airfoil(xs, ys, Re=Re, fileName="tempUAEIII.dat")
#AFtop.Cddata[0, 60:] = 2
#print AFtop.Cddata[0]

"""
print AFroot.Cldata[0]
print AFroot.Cddata[0]
print AFmid.Cldata[0]
print AFmid.Cddata[0]
"""

theTurbine = Turbine([AFroot, AFmid, AFtop], R=R, hubRadius=hubRadius, TSR=TSR, 
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
p2, = pl.plot([3.003, 3.396, 3.760, 4.396, 4.795, 5.168, 5.606, 5.997, 6.394, 6.803, 7.174, 7.594, 7.965, 8.382, 8.788, 9.392, 9.992, 10.482, 10.965, 11.458, 11.980, 12.963, 13.972, 14.975, 15.960, 16.941, 17.946, 18.973, 19.943, 20.957, 21.957], [9.854, 39.740, 89.446, 173.537, 241.123, 311.327, 378.113, 461.949, 549.739, 624.874, 713.140, 823.083, 921.199, 1013.881, 1112.000, 1260.177, 1383.252, 1439.452, 1465.062, 1442.051, 1393.247, 1264.776, 1260.065, 1280.887, 1356.347, 1456.759, 1576.381, 1705.854, 1833.000, 1968.527, 2112.980], 'o')
p3, = pl.plot([2.884, 3.020, 3.156, 3.292, 3.428, 3.564, 3.700, 3.836, 3.972, 4.107, 4.243, 4.379, 4.515, 4.651, 4.787, 4.923, 5.059, 5.195, 5.330, 5.466, 5.602, 5.738, 5.874, 6.010, 6.146, 6.282, 6.418, 6.553, 6.689, 6.825, 6.961, 7.097, 7.233, 7.369, 7.505, 7.641, 7.776, 7.912, 8.048, 8.184, 8.320, 8.456, 8.592, 8.728, 8.864, 9.000, 9.135, 9.271, 9.407, 9.543, 9.679, 9.815, 9.951, 10.087, 10.223, 10.358, 10.494, 10.630, 10.766, 10.902, 11.038, 11.174, 11.310, 11.446, 11.581, 11.717, 11.853, 11.989, 12.125, 12.261, 12.397, 12.533, 12.669, 12.804, 12.940, 13.076, 13.212, 13.348, 13.484, 13.620, 13.756, 13.892, 14.027, 14.163, 14.299, 14.435, 14.571, 14.707, 14.843, 14.979, 15.115, 15.250, 15.386, 15.522, 15.658, 15.794, 15.930, 16.066, 16.202, 16.338, 16.473, 16.609, 16.745, 16.881, 17.017, 17.153, 17.289, 17.425, 17.561, 17.696, 17.832, 17.968, 18.104, 18.240, 18.376, 18.512, 18.648, 18.784, 18.920, 19.055, 19.191, 19.327, 19.463, 19.599, 19.735, 19.871, 20.007, 20.143, 20.278, 20.414, 20.550, 20.686, 20.822, 20.958, 21.094, 21.230, 21.366, 21.501, 21.637, 21.773, 21.909, 22.045, 22.181, 22.317, 22.453], [0, 5.582, 17.756, 29.930, 43.625, 58.842, 75.581, 90.799, 109.059, 128.841, 148.624, 169.928, 192.754, 217.101, 242.970, 270.361, 296.230, 325.143, 354.056, 387.533, 421.011, 454.489, 491.010, 529.053, 565.574, 605.139, 646.225, 688.833, 731.441, 774.049, 816.657, 859.265, 901.873, 944.481, 987.090, 1029.698, 1070.784, 1107.305, 1142.305, 1172.739, 1201.652, 1229.042, 1256.433, 1280.781, 1306.650, 1327.954, 1346.215, 1364.475, 1381.214, 1394.910, 1407.083, 1419.257, 1428.387, 1429.909, 1428.387, 1428.387, 1428.387, 1425.344, 1419.257, 1410.127, 1400.996, 1388.823, 1372.084, 1352.301, 1334.041, 1315.780, 1297.520, 1280.781, 1265.564, 1250.347, 1232.086, 1212.304, 1189.478, 1171.217, 1152.957, 1137.740, 1122.522, 1107.305, 1089.045, 1073.828, 1063.175, 1054.045, 1044.915, 1038.828, 1035.785, 1032.741, 1031.219, 1029.698, 1029.698, 1031.219, 1032.741, 1035.785, 1038.828, 1043.393, 1047.958, 1052.523, 1057.089, 1063.175, 1069.262, 1073.828, 1079.914, 1086.001, 1093.610, 1099.697, 1105.784, 1111.870, 1117.957, 1125.566, 1133.174, 1139.261, 1148.391, 1154.479, 1163.609, 1169.696, 1178.826, 1184.913, 1194.043, 1203.173, 1209.260, 1218.391, 1227.521, 1236.651, 1242.738, 1251.868, 1260.999, 1270.129, 1279.259, 1288.389, 1297.520, 1306.650, 1315.780, 1324.911, 1334.041, 1343.171, 1352.301, 1362.954, 1373.606, 1382.736, 1391.866, 1400.996, 1410.127, 1420.779, 1431.431, 1440.561, 1448.170], '-')

pl.legend([p1, p2, p3], ['Studiens modell (SiteOpt)', 'Experimentell data', u'BEM-data från WT_PERF'], loc=4)

x1,x2,y1,y2 = pl.axis()
pl.axis((0,25,0,2500))

pl.grid()


pl.show()
fig.savefig('Canada.eps', dpi=fig.dpi)


"""
fig4 = pl.figure(figsize=(8, 5))
ax5 = fig4.add_subplot(111)
p1, = pl.plot(AFroot.alfa, AFroot.Cddata[0], linestyle='-')
#p2, = pl.plot(AFS809.alfa, AFS809.Cddata[0], linestyle="--")

x1,x2,y1,y2 = pl.axis()
pl.axis((-3,90,0,2))
#pl.legend([p1], ["asdasdadsdas"], loc=2)
ax5.set_xlabel(r'$\alpha$')
ax5.set_ylabel(u'$C_d$')
pl.grid()
pl.show()

fig5 = pl.figure(figsize=(8, 5))
ax6 = fig5.add_subplot(111)
p1, = pl.plot(AFroot.alfa, AFroot.Cldata[0], linestyle='-')
#p2, = pl.plot(AFS809.alfa, AFS809.Cddata[0], linestyle="--")

x1,x2,y1,y2 = pl.axis()
pl.axis((-3,90,0,1.3))
#pl.legend([p1], ["asdasdadsdas"], loc=2)
ax6.set_xlabel(r'$\alpha$')
ax6.set_ylabel(u'$C_l$')
pl.grid()
pl.show()
"""

TSRs = [2,3,4,5,6,7,8,9,10,11,12,13] 
Re = np.array([1])

Cps = []

cutIn = 3.6
cutOut = 3.6
RPM = None

for TSR in TSRs:

	theTurbine = Turbine([AFroot, AFmid, AFtop], R=R, hubRadius=hubRadius, TSR=TSR, 
	                      ratedPower=ratedPower, B=B, visc=visc, rho=rho, tol=tol, 
	                      windSpeeds=windSpeeds, windFrq=windFrq, cutIn=cutIn, 
	                      cutOut=cutOut, skipBlending=1, RPM=RPM,
	                      chord=chord, twist=twist, pitch=pitch)
	Cps.append(theTurbine.Cps[0])
	print theTurbine.Cps[0]


fig4 = pl.figure(figsize=(8, 5))
ax4 = fig4.add_subplot(111)
p1, = pl.plot(TSRs, Cps)
p2, = pl.plot([1.517, 1.670, 1.847, 2.093, 2.216, 2.374, 2.568, 2.776, 2.905, 3.031, 3.185, 3.343, 3.494, 3.633, 3.896, 4.085, 4.284, 4.522, 4.769, 5.062, 5.382, 5.561, 5.762, 6.210, 6.423, 6.700, 6.981, 7.280, 7.980, 8.348, 8.798, 9.309, 10.462, 11.150], [0.037, 0.044, 0.052, 0.063, 0.073, 0.088, 0.110, 0.154, 0.182, 0.211, 0.237, 0.263, 0.279, 0.296, 0.319, 0.336, 0.351, 0.364, 0.372, 0.388, 0.408, 0.407, 0.409, 0.419, 0.424, 0.423, 0.416, 0.403, 0.370, 0.352, 0.320, 0.264, 0.136, 0.023], 'o')

x1,x2,y1,y2 = pl.axis()
#pl.axis((-3,25,0,0.25))
pl.legend([p1, p2], ["Studiens modell (SiteOpt) 3,6 m/s", "Experimentell data 3,6 m/s"], loc=8)
ax4.set_xlabel(r'Tip Speed Ratio $\lambda$')
ax4.set_ylabel(u'$C_p$')
pl.grid()
pl.show()
