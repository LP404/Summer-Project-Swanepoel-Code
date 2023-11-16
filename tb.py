import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import csv
from scipy.optimize import curve_fit


x = np.array([2,2,2])
d = np.array([2,2,2])

def ycalc(d):
    for i in range(len(d)):
        k = x + d
    
    return k


def Richards_Penman(xIn,A,K,B,v,Q,M,C,L,P):
    final = A + (K-A+(L*np.log((xIn-M+P))))/((C + np.exp(-B*(xIn-M)))**(v))
    return final

def logifunc(x,A,x0,k,off):
    return A / (1 + np.exp(-k*(x-x0)))+off


def FancyInterpolate(xInterp,xArray,yArray,AntiNodeX,AntiNodeY,BoolMinima):
    
    #Loc0 = np.where((np.diff(yArray) / max(np.diff(yArray))) == 1)[0][0]
    
    if BoolMinima == True:
        Loc = max(np.where(yArray < AntiNodeY[0])[0])
    else: 
        Loc = np.where(xArray ==  AntiNodeX[0])[0][0]
        
    #Loc = max(np.where((np.around(yArray,2) == 0.50) | (np.around(yArray,2) == 0.49) | (np.around(yArray,2) == 0.51))[0])
    
    #0.5 is currently arbiraty, but frindges don't usually appear until after 50% transmission
    
    x1 = np.append(xArray[0:Loc],AntiNodeX)
    y1 = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(x1)
    x1 = x1[Sort[::1]]
    y1 = y1[Sort[::1]]
        
    yInterp = np.interp(xInterp,x1,y1)
    
    return yInterp

def FancyInterpolate2(xInterp,xArray,yArray,AntiNodeX,AntiNodeY,BoolMinima):
    
    #Loc0 = np.where((np.diff(yArray) / max(np.diff(yArray))) == 1)[0][0]
    
    # if BoolMinima == True:
    #     Loc = max(np.where(yArray < AntiNodeY[0])[0])
    # else: 
    #     Loc = np.where(xArray ==  AntiNodeX[0])[0][0]
        
    Loc = max(np.where((np.around(yArray,2) == 0.60) | (np.around(yArray,2) == 0.59) | (np.around(yArray,2) == 0.61))[0])
    
    #0.5 is currently arbiraty, but frindges don't usually appear until after 50% transmission
    
    x1 = np.append(xArray[0:Loc],AntiNodeX)
    y1 = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(x1)
    x1 = x1[Sort[::1]]
    y1 = y1[Sort[::1]]
        
    guess = [0,0.5,50,10,1,200,1,0.25,1]
    
    yInterp = np.interp(xInterp,x1,y1)
    
    
    
    TransFit, Resi = curve_fit(Richards_Penman,x1,y1,guess, maxfev=500000000)
    yInterp2 = Richards_Penman(xInterp,TransFit[0],TransFit[1],TransFit[2],TransFit[3],TransFit[4],TransFit[5],TransFit[6],TransFit[7],TransFit[8])
    
    return yInterp, yInterp2, x1,y1


#dt = 1

h = 6.63e-34
c = 3e8

#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))
# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData\\LabData\\DP'))

#In nm
Slitwidth = 2

xMin = 200
xMax = 800

EnableCorrection = True

#Disable multiplier before using
EnableCorrection = True

A1 = 1.4313493
A2 = 0.65054713
A3 = 5.3414021
B1 = 5.2799261e-3
B2 = 1.42382647e-2
B3 = 325.017834

#Transmission unceranity
TU = (0.5/100)
#In nm
DeltaLam = 0.3



#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
suffix = ".txt"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)


    #Change 0.075 back to 0.2
    vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])
    
    
    if i == 0:
        #All the xValues are the same, only need to do this once
        xP = np.linspace(min(vars()[files[i]+'T'][0]),max(vars()[files[i]+'T'][0]),10001)
    else:
        pass

for i in range(len(files1)):
    files1[i] = files1[i].rstrip(".txt")
    
    vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    
    #This can be commented out with an alternet approach
    vars()[files1[i]][0] = vars()[files1[i]][0]
    
    vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
    # vars()[files1[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files1[i]],xMin,xMax))
    vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],xMin,xMax)

    vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])

    vars()['nS'+files1[i]] = F.Sub_nFinder(vars()['yP'+files1[i]])

if EnableCorrection == True:    
    for i in range(len(files)):
        
        MaxLocOriginal = np.where(vars()['yFiltered'+str(i)] == max(vars()['yFiltered'+str(i)]))[0][0]
        
        Val = vars()[files[0]+'T'][0][MaxLocOriginal]
        
        MaxLocSub = np.where(xP == F.FindNearestVal(xP,Val))[0][0]
        
        SubValue = vars()['yP'+files1[0]][MaxLocSub]
        
        CorrectiveMultiplier = vars()['yP'+files1[0]][MaxLocSub] / vars()['yFiltered'+str(i)][MaxLocOriginal]
        
        vars()['yFiltered'+str(i)] =  vars()['yFiltered'+str(i)] * CorrectiveMultiplier

else:
    pass    

for i in range(len(files)):            
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])
    
    
    #0.075 is standard value for the second variable in the function. but adjust as neccessary
    vars()['xNewMax'+str(i)], vars()['yNewMaxUnCorr'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMinUnCorr'+str(i)] = F.BackNodeFix(3,0.075,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)])
    
    vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)]
    vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)]

    vars()['yPMaxUnCorr'+str(i)], vars()['yPMaxUnCorrRfit'+str(i)], vars()['yPMaxxRawUncor'+str(i)], vars()['yPMaxyRawUncor'+str(i)]  = FancyInterpolate2(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMaxUnCorr'+str(i)],False)
    vars()['yPMinUnCorr'+str(i)], vars()['yPMinUnCorrRfit'+str(i)], vars()['yPMinxRawUncor'+str(i)], vars()['yPMinyRawUncor'+str(i)]  = FancyInterpolate2(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMinUnCorr'+str(i)],True)

    vars()['yPMax'+str(i)], vars()['yPMaxRfit'+str(i)], vars()['yPMaxxRaw'+str(i)], vars()['yPMaxyRaw'+str(i)] = FancyInterpolate2(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)],False)
    vars()['yPMin'+str(i)], vars()['yPMinRfit'+str(i)], vars()['yPMinxRaw'+str(i)], vars()['yPMinyRaw'+str(i)]  = FancyInterpolate2(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)],True)
    
    np.savetxt(files[i]+'xMaxUncor.txt',vars()['yPMaxxRawUncor'+str(i)])
    np.savetxt(files[i]+'yMaxUncor.txt',vars()['yPMaxyRawUncor'+str(i)])
    np.savetxt(files[i]+'xMinUncor.txt',vars()['yPMinxRawUncor'+str(i)])
    np.savetxt(files[i]+'yMinUncor.txt',vars()['yPMinyRawUncor'+str(i)])
    np.savetxt(files[i]+'xMax.txt',vars()['yPMaxxRaw'+str(i)])
    np.savetxt(files[i]+'yMax.txt',vars()['yPMaxyRaw'+str(i)])
    np.savetxt(files[i]+'xMin.txt',vars()['yPMinxRaw'+str(i)])
    np.savetxt(files[i]+'yMin.txt',vars()['yPMinyRaw'+str(i)])


for i in range(len(files)):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.figure(i,figsize = (8,5),dpi = 600)
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)], label = "Filtered Data")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    
    # np.savetxt(str(files[i])+".csv", vars()[files[i]+'T'].T, delimiter=',')


    plt.scatter(vars()['xNewMax'+str(i)] , vars()['yNewMaxUnCorr'+str(i)], color = 'black', marker = "x", label = "Maxima")
    plt.scatter(vars()['xNewMin'+str(i)] , vars()['yNewMinUnCorr'+str(i)], color = 'red', marker = "x", label = "Minima")
    plt.plot(xP,vars()['yPMaxUnCorr'+str(i)], color = 'black', linestyle="dotted", label = "TM")
    plt.plot(xP,vars()['yPMinUnCorr'+str(i)], color = 'red', linestyle="dotted", label = "Tm")
    plt.plot(xP,vars()['yPMaxUnCorrRfit'+str(i)], color = 'blue', linestyle="dashed", label = "TM Rfit")
    plt.plot(xP,vars()['yPMinUnCorrRfit'+str(i)], color = 'green', linestyle="dashed", label = "Tm Rfit")
     
    
    plt.plot(xP,vars()['yP'+files1[0]], color = 'orange', label = "Substrate")
    
    plt.legend()
    
