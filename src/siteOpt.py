#!/usr/bin/env python
# encoding: utf-8
# python v 2.7
"""
siteOpt.py

Created by Max Berggren 2014-05-10

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

from __future__ import division
import numpy as np
from scipy import interpolate
import subprocess as sp
from StringIO import StringIO
from math import pi, atan, acos, cos, sin, tan, exp, sqrt, radians, degrees
import os

xfoilPath = os.getcwd() + os.path.sep + "Xfoil.app/Contents/Resources/xfoil"


def BlendAirfoils(AF1, AF2, percentAF2):
    delX = AF2.x - AF1.x
    newx = AF1.x + delX*percentAF2
    delY = AF2.y - AF1.y
    newy = AF1.y + delY*percentAF2
    return newx, newy

def AoAsteps(start, stop, step):
    AoAs = []
    AoA = start
    if start <= stop:
        while AoA <= stop:
            AoAs.append(AoA)
            AoA += step
    else:
        while AoA >= stop:
            AoAs.append(AoA)
            AoA -= step
    return AoAs

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm


def individualNr():
    counter = 0  
    while True:
        counter = counter + 1
        yield counter

def readTUR(fileName):
    with open(fileName) as f:
        content = f.readlines()
        content = "   ".join(content)
        content = content.replace("\n", "")
        content = content.replace("[", " ")
        
        content = content.replace("]", " ")
        content = content.replace("     ", " ")
        content = content.replace("    ", " ")
        content = content.replace("   ", " ")
        content = content.replace("  ", " ")
        content = content.split(" ")
        
        content = [x for x in content if x]
        content = [float(x) for x in content]
        content = np.array(content)
       
        
        return content


class Turbine:
    def __init__(self, listOfAFs, R, hubRadius, TSR, ratedPower, B, visc, rho, tol, windSpeeds, 
                 windFrq, cutIn=5, cutOut=15, noBetween=2, noInterpolatedAFs=1, 
                 skipBlending=1, metaData=None, indNr=None, RPM=None,
                 chord=[0.4572,  0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572, 0.4572], 
                 twist=[44, 44, 36, 24, 20, 16, 13, 10, 9, 4, 2, 1, 0], pitch=0,
                 chordData=None, twistData=None,
                 evalWindSpeeds=None):

        self.radialPositions = len(listOfAFs) + (len(listOfAFs)-1)*noBetween + ((len(listOfAFs)-1)*noBetween-1+len(listOfAFs))*noInterpolatedAFs

        self.AFs = [None] * self.radialPositions # initiate list of AFs
        self.R = R
        self.hubRadius = hubRadius # radius where the blade begins and hub ends (m)
        self.TSR = TSR # tip speed ratio (a.k.a. the symbol lambda)
        self.RPM = RPM
        self.B = B # number of blades
        self.visc = visc # kinematic viscosity of air 
        self.rho = rho # density of air 
        self.tol = tol # convergence tolerance for a and aprime
        self.cutIn = cutIn
        self.cutOut = cutOut
        self.windSpeeds = windSpeeds
        self.windFrq = windFrq
        self.ratedPower = ratedPower
        self.chord = chord
        self.twist = twist
        self.pitch = pitch
        self.r = []

        if evalWindSpeeds == None:
            self.evalWindSpeeds = np.linspace(cutIn, cutOut, 10)
        else:
            self.evalWindSpeeds = evalWindSpeeds

        if chordData:
            self.chord = self.interpolChordTwist(chordData)
        if twistData:
            self.twist = self.interpolChordTwist(twistData)
        if indNr:
            self.writeMetaData(metaData, indNr)

        j = 0
        for i in range(len(listOfAFs)): # portion out the main AFs
            self.AFs[int(i*((self.radialPositions-1)/(len(listOfAFs)-1)))] = listOfAFs[j]
            j = j+1


        for i in range(len(listOfAFs)-1): # find the places of the blended AFs
            startIndex = int(i*((self.radialPositions-1)/(len(listOfAFs)-1)))
            nextMainAF = int((i+1)*((self.radialPositions-1)/(len(listOfAFs)-1)))
            betweenStep = int((self.radialPositions/len(listOfAFs))/(noBetween))
            h = 1
            for b in range(noBetween+1):
                
                if not b == 0:
                    if not skipBlending:
                        xb, yb = BlendAirfoils(listOfAFs[i], listOfAFs[i+1], percentAF2=h*0.3333333)
                        AFb = Airfoil(xb, yb, blendInput=1)
                        self.AFs[int(startIndex + b*betweenStep)] = AFb
                    else:
                        self.AFs[int(startIndex + b*betweenStep)] = Airfoil(interpolatedInput=1,
                                                                            AF1=listOfAFs[i], 
                                                                            AF2=listOfAFs[i+1])

        for i in range(len(listOfAFs)-1): # find the places of the intepolated AFs,
            startIndex = int(i*((self.radialPositions-1)/(len(listOfAFs)-1)))
            betweenStep = int((self.radialPositions/len(listOfAFs))/(noBetween))
            for b in range(noBetween+1):
                for k in range(noInterpolatedAFs):
                    #print "interp ",startIndex + b*betweenStep + k +1
                    interIndex = int(startIndex + b*betweenStep + k +1)
                    self.AFs[interIndex] = Airfoil(interpolatedInput=1, AF1=self.AFs[interIndex-1], AF2=self.AFs[interIndex+1])


        self.checkPerformance() # run simulation

    def interpolChordTwist(self, data):
        tck = interpolate.splrep([0, data[1], 1], [data[0], data[2], data[3]],s=0,k=2)
        xnew = np.linspace(0, 1, self.radialPositions)
        ynew = interpolate.splev(xnew,tck,der=0)        
        return ynew
    
    def writeMetaData(self, metaData, indNr):
        with open("TurbineNo"+str(indNr)+".tur", "w+") as f:
            f.write(str(metaData) + '\n') 

    def avgPower(self):
        if self.power:
            avgPower = 0

            for i, speed in enumerate(self.windSpeeds):
                if speed > self.cutIn and speed < self.cutOut:
                    avgPower += self.windFrq[i]*self.powerInterpolate(speed)
            return avgPower

    def powerInterpolate(self, windSpeed):
        try:
            if windSpeed < self.cutIn:
                return 0
            if windSpeed > self.cutOut:
                return 0
            
            if len(self.evalWindSpeeds) == 1:
                return self.power[0]  

            f = interpolate.interp1d(self.evalWindSpeeds, self.power, kind='cubic')
            if f(windSpeed) < self.ratedPower:
                return f(windSpeed)
            else:
                return self.ratedPower
        except:
            print "Problems with power curve"
            return 0

    def checkPerformance(self):
        R = self.R # radius (m)
        hubRadius = self.hubRadius # radius where the blade begins and hub ends (m)
        TSR = self.TSR # tip speed ratio (a.k.a. the symbol lambda)
        B = self.B # number of blades
        visc = self.visc # kinematic viscosity of air 
        rho = self.rho # density of air 
        tol = self.tol # convergence tolerance for a and aprime
        cutIn = self.cutIn
        cutOut = self.cutOut
        RPM = self.RPM
        pitch = self.pitch

        r = np.linspace(hubRadius/R, 1, num=self.radialPositions)*R # radius of blade elements betweeb hub and tip (m)
        self.r = r
        dr = r[1]-r[0] # width of blade elements

        rNorm = r/R # normalized radial elements (r/R)

        twist = np.array(self.twist) + pitch # twist distribution accounted for pitch
        chord = np.array(self.chord) # chord distribution

        T = 0.0
        Q = 0.0
        Cp = 0

        powers = []
        Qs = []
        angSpeeds = []
        Cps = []

        

        for Uinf in self.evalWindSpeeds:
            
            angSpeed = TSR*Uinf/R # radians/s
            if RPM: # if constant RPM is to be used
                angSpeed = RPM*0.1047 # radians/s
                TSR = angSpeed*R/Uinf

            angSpeeds.append(angSpeed*(9.5493)) # rpm
            T = 0
            Q = 0

            print "Uinf    Element Radius  Iter    a       aprime  AoA     Utot    locTSR  phi     Cl      Cd      F       sigma   Re*10^5 AngSpd  phi(deg)  dQ    Cn      Ct"
            for i in range(len(r)-1): 

                a = 0
                aprime = 0
                dQ = 0
                dT = 0
                
                TSRr = TSR*rNorm[i] # local TSR at the blade element (a.k.a. lambda_r)
                sigmaprime = B*chord[i]/(2*pi*r[i]) # local solidity


                #a = (1/4)*(2 + pi*TSRr*sigmaprime - sqrt(4 - 4*pi*TSRr*sigmaprime + pi*sigmaprime**2*(8*twist[i] + pi*sigmaprime)))


                for j in range(200): # max 200 iterations
                
                    if not aprime == -1:
                        phi = atan( (1 - a)/( (1 + aprime)*TSRr ) )
               
                    else: # prevent divide by zero
                        phi = atan( (1 - a)/( (1 + aprime*1.01)*TSRr ) )

                    AoA = phi*(180.0/pi) - twist[i] # degrees
                    Utot = sqrt(  (Uinf*(1 - a))**2 + ( angSpeed*r[i]*(1 + aprime) )**2  )
                    Re = Utot*chord[i]/visc
                    Cl = self.AFs[i].Cl(AoA, Re/100000)
                    Cd = self.AFs[i].Cd(AoA, Re/100000)


                    Cn = Cl*cos(phi) + Cd*sin(phi)
                    Ct = Cl*sin(phi) - Cd*cos(phi)

                    
                    C_T = sigmaprime*(1 - a)**2*Cn/(sin(phi)**2)

                    try:
                        Ftip = (2/pi)*acos(exp(-(B*(R-r[i]))/(2*r[i]*sin(phi))))
                        Rhub = r[0]
                        Fhub = (2/pi)*acos(exp(-(B*(r[i]-Rhub))/(2*r[i]*sin(phi))))
                        F = Fhub*Ftip
                    except: # First element often fails
                        F = 0.5
                    
                    np.seterr(all='raise') 
                    
                    if C_T > 0.96*F: # Glauert correction
                        newa = (18*F - 20 - 3*sqrt(C_T*(50-36*F) + 12*F*(3*F - 4))) / (36*F - 50)
                    else: # Standard BEM therory
                        try:
                            newa = 1 / ( 1 + (4*F*sin(phi)**2)/(sigmaprime*Cn) )
                        except: # prevent divide by zero
                            print "error när a räknas ut"
                            dQ = 0
                            dT = 0
                            newa = 0
                            break
                    try:
                        aprime = 1 / (-1 + (4*F*sin(phi)*cos(phi)) / (sigmaprime*Ct))
                    except: # prevent divide by zero
                        dQ = 0
                        dT = 0
                        aprime = 0
                        break


                    diffa = abs(a - newa)
                    damper = 0.5
                    a = damper*newa + (1 - damper)*a

                    if diffa < tol and j > 3:
                        break
                
                #dQ = 4*pi*(r[i]**3)*rho*Uinf*angSpeed*(1 - a)*aprime*dr
                #dQ = 0.5*rho*B*Uinf**2*(1-a)*angSpeed*r[i]*(1+aprime)*chord[i]*(Cl*sin(phi) - Cd*cos(phi)*r[i]*dr/(sin(phi)*cos(phi))
                #dQ = B*0.5*rho*(Utot**2)*Ct*chord[i]*r[i]*dr
                #dQ = 4*aprime*(1 - a)*rho*Uinf*angSpeed*r[i]**3*pi*dr

                dQ = 0.5*rho*B*Uinf*(1-a)*angSpeed*r[i]*(1+aprime)*chord[i]*Ct*r[i]*dr/(sin(phi)*cos(phi)) # calc additional Q
                
                if j == 199: # convergence fail
                    dQ = 0
                    dT = 0

                if AoA < -20 or AoA > 90: # if crazy AoA
                    dQ = 0
                    dT = 0

                Q = Q + dQ # total torQue

                print str(int(Uinf))+"\t"+str(i)+"\t"+str('%.2f' % r[i])+"\t"+str(j)+"\t"+str('%.2f' % a)+"\t"+str('%.2f' % aprime)+"\t"+str('%.2f' % AoA)+"\t"+str('%.2f' % Utot)+"\t"+str('%.2f' % TSRr)+"\t"+str('%.2f' % phi)+"\t"+str('%.2f' % Cl)+"\t"+str('%.2f' % Cd)+"\t"+str('%.2f' % F)+"\t"+str('%.2f' % sigmaprime)+"\t"+str('%.2f' % (Re/100000))+"\t"+str('%.2f' % angSpeed)+"\t"+str('%.2f' % (phi*180/pi))+"\t"+str('%.2f' % dQ)+"\t"+str('%.2f' % Cn)+"\t"+str('%.2f' % Ct)



            Power = Q*angSpeed
            Cp = Power / (0.5*rho*Uinf**3*pi*R**2)
            powers.append(Power)
            Qs.append(Q)
            Cps.append(Cp)


        self.power = powers
        self.torque = Qs
        self.RPMs = angSpeeds
        self.Cps = Cps




class Airfoil:

    def __init__(self, x=None, y=None, Re=[0.5, 1, 3, 4, 5, 7, 10, 20, 50, 70, 100, 140], 
                             AoAstart=-3, 
                             AoAstop=25, 
                             AoAstep=2, 
                             Ncrit=9,
                             blendInput=0,
                             interpolatedInput=0,
                             AF1 = None,
                             AF2 = None,
                             fileName="temp.dat",
                             highAngleInterp=1,
                             alfa=None,
                             Cldata=None,
                             Cddata=None,
                             AR=None):

        if not interpolatedInput:
            self.interpolatedInput = False
            self.originalx = x
            self.originaly = y

            if not blendInput:
                self.x, self.y = np.append(x, 1), np.append(y, 0) # add the constant front
                self.x, self.y = np.insert(self.x, 0, 1), np.insert(self.y, 0, 0) # and back coordinates
                self.x, self.y = self.connectTheDots(self.x, self.y) # interpolate between
            else:
                self.x, self.y = x, y


            self.alfa = AoAsteps(AoAstart, AoAstop, AoAstep)
            self.Cldata = []
            self.Cddata = []
            self.Cmdata = []
            self.Re = []
            self.AoAstop = AoAstop
            self.AoAstart = AoAstart
            self.AR = AR

            if Cldata and Cddata: # if empirical data is to be used instead of XFOIL
                if not len(Re) == 1:
                    print "If empirical Cl- and Cddata is used, only one Reynolds can be used matching the data."
                else:
                    self.Cldata.append(Cldata)
                    self.Cddata.append(Cddata)
                    self.alfa = alfa
                    self.AoAstop = max(self.alfa)
                    self.AoAstart = min(self.alfa)
                    self.Re = Re

            else: # if xfoil is to be used

                self.writeDAT(self.x, self.y, fileName) # write to disk so XFOIL can read them
                
                for i, Rei in enumerate(Re): 
                    
                    alfa, Cl, Cd, Cm = self.getPolars(Rei*10**5, 
                                                      AoAstart, 
                                                      AoAstop, 
                                                      AoAstep, 
                                                      Ncrit,
                                                      airfoil=fileName)
                    
                    if Cl:
                        self.Cldata.append(Cl)
                        self.Cddata.append(Cd)
                        self.Cmdata.append(Cm)
                        self.Re.append(Rei)

            try:
                self.Cldata = np.array(self.Cldata)
                self.Cddata = np.array(self.Cddata)
                self.Cmdata = np.array(self.Cmdata)

                self.Cldata = self.fillHoles(self.Cldata)
                self.Cddata = self.fillHoles(self.Cddata)
                self.Cmdata = self.fillHoles(self.Cmdata)

                if highAngleInterp:
                    self.highAngleInterp()
                    
                # OBS: kass fix kika närmare på!
                #self.Cddata = self.Cddata*2

            except:
                self.failed = True
            
        else: # interpolate between two existing AFs Cl and Cd
            self.interpolatedInput = True
            self.AF1 = AF1
            self.AF2 = AF2
            self.alfa = AF1.alfa
            try:
                self.AR = (AF1.AR + AF2.AR)/2
            except:
                if AF1.AR:
                    self.AR = AF1.AR
                else:
                    self.AR = AF2.AR


    def highAngleInterp(self):
        ClnewdataR = []
        CdnewdataR = []
        ClnewdataL = []
        CdnewdataL = []
        CLtemp = []
        CDtemp = []

        for i, each in enumerate(self.Cldata):
            # Todo: fix!
            R = 10.046/2
            chord = 0.4572
            AR = R/chord
            if self.AR:
                AR = self.AR
            print "AR=",AR
            if AR >= 50:
                Cdmax = 2
            else:
                Cdmax = 1.11 + 0.018*AR
            B1 = Cdmax

            stallAlpha = radians(self.alfa[np.argmax(self.Cldata[i])])   
            stallAlphaInd = np.argmax(self.Cldata[i])    
            Clstall = self.Cldata[i,stallAlphaInd]
            Cdstall = self.Cddata[i,stallAlphaInd]
            B2 = (Cdstall - Cdmax*sin(stallAlpha)**2)/cos(stallAlpha)

            A1 = B1/2
            A2 = (Clstall - Cdmax*sin(stallAlpha)*cos(stallAlpha))*sin(stallAlpha)/(cos(stallAlpha)**2)

            # from alfa-stall to end of vector
            alphas = np.radians(self.alfa[stallAlphaInd:])
            alphas[alphas == 0] = 0.1 # to prevent from divide by zero

            CD = B1*np.sin(alphas)**2+B2*np.cos(alphas)
            CL = A1*np.sin(2*alphas) + A2*np.cos(alphas)**2/np.sin(alphas)

            self.Cldata[i,stallAlphaInd:] = CL
            self.Cddata[i,stallAlphaInd:] = CD

            # from end of vector to alfa = 90
            anglesTo90 = np.linspace(self.AoAstop+1, 89, 40)
            anglesTo90rad = np.radians(anglesTo90)

            CD = B1*np.sin(anglesTo90rad)**2 + B2*np.cos(anglesTo90rad)
            CL = A1*np.sin(2*anglesTo90rad) + A2*np.cos(anglesTo90rad)**2/np.sin(anglesTo90rad)
            
            # alfa = 90 to 180
            angles90to180 = np.linspace(90, 180, 40)
            angles90to180rad = np.radians(angles90to180)

            CL = np.hstack((CL, 2*np.sin(angles90to180rad)*np.cos(angles90to180rad)))
            CD = np.hstack((CD, B1*np.sin(angles90to180rad)**2))
            #CD = np.hstack((CD, 2*np.sin(angles90to180rad)**2)) WHUT?
            
            # from alfa = -180 to start of vector 
            anglesNeg180toNeg90 = np.linspace(-180, -45, 40)
            anglesNeg180toNeg90rad = np.radians(anglesNeg180toNeg90)
            
            CLneg18090 = 2*np.sin(anglesNeg180toNeg90rad)*np.cos(anglesNeg180toNeg90rad)
            CDneg18090 = B1*np.sin(anglesNeg180toNeg90rad)**2
            
            ClnewdataR.append(CL)
            CdnewdataR.append(CD)

            ClnewdataL.append(CLneg18090)
            CdnewdataL.append(CDneg18090)

        self.Cldata = np.hstack((ClnewdataL, self.Cldata, np.array(ClnewdataR)))
        self.Cddata = np.hstack((CdnewdataL, self.Cddata, np.array(CdnewdataR)))

        self.alfa = np.hstack((anglesNeg180toNeg90, self.alfa, anglesTo90, angles90to180))


    def fillHoles(self, data):
        """ filling holes where XFOIL didn't converge instead of running again """
        for i, each in enumerate(data):

            try:
                firstInd = np.nonzero(each)[0][0]
                lastInd = np.nonzero(each)[0][-1]

                fix = each[firstInd:lastInd]
                x, = np.nonzero(fix)

                fixed = np.interp(np.arange(len(fix)), x, fix[x])

                data[i, firstInd:lastInd] = fixed
            except:
                """
                Probably just zeros
                """

        return data

    def connectTheDots(self, x, y):
        tck,u = interpolate.splprep([x,y],s=0) # Find the apropiate spline
        unew = np.arange(0,1.005,0.005)
        out = interpolate.splev(unew,tck)

        return out[0], out[1]
    
    def Cl(self, AoA, Re):
        if self.interpolatedInput: # return halfway of two other AFs
            try:
                return (self.AF1.Cl(AoA, Re) + self.AF2.Cl(AoA, Re)) / 2
            except:
                return 0
        else:
            try:
                if len(self.Re) == 1:
                    f = interpolate.interp1d(self.alfa, self.Cldata[0], kind='linear')
                    return f(AoA)
                else:
                    f = interpolate.interp2d(self.alfa, self.Re, self.Cldata, kind='linear')
                    return f(AoA, Re)[0]
                
            except:
                return 0

    def Cd(self, AoA, Re):
        if self.interpolatedInput: # return halfway of two other AFs
            try:
                return (self.AF1.Cd(AoA, Re) + self.AF2.Cd(AoA, Re)) / 2
            except:
                return 2
        else:
            try:
                if len(self.Re) == 1:
                    f = interpolate.interp1d(self.alfa, self.Cddata[0], kind='linear')
                    return f(AoA)
                else:
                    f = interpolate.interp2d(self.alfa, self.Re, self.Cddata, kind='linear')
                    return f(AoA, Re)[0]
            except:
                return 2

    def writeDAT(self, x, y, fileName="temp.dat"):
        with open(fileName, "w+") as f:
            decimals = 4
            f.write(fileName + '\n') # Write fileName at very top of file
            for i in range(len(x)):
                f.write(str(round(x[i],decimals)) + "     " + str(round(y[i],decimals)) + '\n')

    def getPolars(self, Re, AoAstart, AoAstop, AoAstep, Ncrit=9, airfoil="temp.dat", surpressGUI=1):

        def issueCmd(cmd, echo=True):
            ps.stdin.write(cmd + '\n')

        ps = sp.Popen([xfoilPath], 
            stdin=sp.PIPE, 
            stdout=sp.PIPE,
            stderr=None,
            shell=True)

        try:
            os.remove(str(airfoil) + '.pol') # remove file if it already exists
        except:
            """
            #print "Kunde inte hitta gammal polare"
            """
        issueCmd('load ' + str(airfoil))
 
        if surpressGUI: # make XFOIl surpress the visuals since they can make you go cray cray
            issueCmd('PLOP')
            issueCmd('G')
            issueCmd('')

        issueCmd('PANE') # adds points if needed
        issueCmd('PANE')
        issueCmd('OPER')

        if not Ncrit == 9:
            issueCmd('vpar')
            issueCmd('n ' + str(Ncrit))
            issueCmd('')

        issueCmd('VISC ' + str(Re))

        issueCmd('iter 50')
        issueCmd('PACC')
        issueCmd(str(airfoil) + '.pol')
        issueCmd('')
        issueCmd('ASEQ '+str(AoAstart)+' '+str(AoAstop)+' '+str(AoAstep))
        #issueCmd('ASEQ -2.5 -2.0 0.05')
        #issueCmd('ASEQ -1.5  8.0 0.5')
        #issueCmd('ASEQ  8.2  9.0 0.2')
        issueCmd('PACC')
        issueCmd('')
        issueCmd('quit')
        outputFromTerminal = ps.stdout.read()
        #print outputFromTerminal

        with open(str(airfoil) + '.pol') as f: # read file from XFOIL
            try:
                content = f.readlines()
                if len(content) > 14:
                    content = StringIO("\n".join(content[12:]))
                    content = np.loadtxt(content)
                    alfa = content[:,0]
                    CL = content[:,1]
                    CD = content[:,2]
                    CM = content[:,3]

                    CLDict = dict(zip(alfa, CL))
                    CDDict = dict(zip(alfa, CD))
                    CMDict = dict(zip(alfa, CM))

                    CL, CD, CM = [], [], []

                    for AoA in AoAsteps(AoAstart, AoAstop, AoAstep):
                        try:
                            CL.append(CLDict[AoA])
                            CD.append(CDDict[AoA])
                            CM.append(CMDict[AoA])
                        except KeyError:
                            CL.append(0)
                            CD.append(0)
                            CM.append(0)

                    print alfa, CL, CD, CM
                    return alfa, CL, CD, CM
                else:
                    return None, None, None, None
                 
            except:
                return None, None, None, None            


def readDAT(fileName):
    with open(fileName) as f:
        x, y = [], []
        content = f.readlines()

        for row in content:
            try:
                left = float(row.split("     ")[0])
                right = float(row.split("     ")[1].replace("\r\n",""))
                
                x.append(left)
                y.append(right)
            except:
                """
                Probably a name of the airfoil at top of file
                """
        return np.array(x), np.array(y)  

