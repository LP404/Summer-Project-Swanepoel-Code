import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import Functions1 as F1
import csv
from scipy.optimize import curve_fit

def Reinterp(xInterp,xArray,yArray,AntiNodeX,AntiNodeY,CutOff):
    #Since we know the cutoff is 0.6
    #Given we know how the array is constructed we can assume the next index is for the first maxima/minima

    Loc = max(np.where((np.around(yArray,2) == CutOff) | (np.around(yArray,2) == CutOff-0.01) | (np.around(yArray,2) == CutOff+0.01))[0])
    
    xNuvo = np.append(xArray[0:Loc],AntiNodeX)
    yNuvo = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(xNuvo)
    xNuvo = xNuvo[Sort[::1]]
    yNuvo = yNuvo[Sort[::1]]

    yNuvoInd = np.where(yNuvo == F.FindNearestVal(yNuvo,CutOff))[0][0]
    y1 = yNuvo[yNuvoInd - 4]
    x1 = xNuvo[yNuvoInd - 4]
    y2 = yNuvo[yNuvoInd + 1]
    x2 = xNuvo[yNuvoInd + 1]
    x0, y0, a, b, c = F.ellipseParamFinder(x1,y1,x2,y2)   
    
    xEllipse = np.arange(x1, x2+1, 1)
    yEllipse = F.quarterEllipseFinder(xEllipse,a,b,x0,y0)
        
    guess = np.array([np.exp(1),1,260])
    
    yTopStart =   yNuvo[yNuvoInd+1:]
    xTopStart =   xNuvo[yNuvoInd+1:]
      
    TransFit, Resi = curve_fit(F.logCurve,xTopStart, yTopStart,guess, maxfev=500000000)
       
    xTopStart2 = np.arange(xNuvo[yNuvoInd +1], 801, 1)
        
    yFitbutIDK = F.logCurve(xTopStart2 , TransFit[0], TransFit[1], TransFit[2])

    #Since step size is 1 we can use np.gradient
    
    #!!! Note: Change X values to cover all x points for final fit
    
    gradEllipse = np.gradient(yEllipse)
    bEllipse = yEllipse - (gradEllipse*xEllipse)
    
    gradyFitbutIDK = np.gradient(yFitbutIDK)
    bFitbutIDK = yFitbutIDK - (gradyFitbutIDK*xTopStart2)
 
    allowance = 20

    xArr = np.linspace(x0-allowance,x0+allowance,100001)
    IncLoc = F.LineComp(xArr,gradEllipse,gradyFitbutIDK,bEllipse,bFitbutIDK,5)


    yArr = (xArr*gradEllipse[IncLoc[1]] + bEllipse[IncLoc[1]])
    yArr1 = (xArr*gradyFitbutIDK[IncLoc[2]] + bFitbutIDK[IncLoc[2]])
    
    xArr = np.round(xArr).astype(int)

    GeoMeanNewyArr = np.sqrt((yArr * yArr1))
    
    newxArr, newGeoMeanNewyArr = F.ArrayCondenser(xArr,GeoMeanNewyArr)

    np.savetxt('ErrorLogx0.txt',np.array([x0]))
    np.savetxt('ErrorLogxEllipse.txt',xEllipse)
    np.savetxt('ErrorLogyEllipse.txt',yEllipse)
    np.savetxt('ErrorLogxTopStart2.txt',xTopStart2)
    np.savetxt('ErrorLogyFitbutIDK.txt',yFitbutIDK)
    np.savetxt('ErrorLognewxArr.txt',newxArr)
    np.savetxt('ErrorLognewGeoMeanNewyArr.txt',newGeoMeanNewyArr)
    np.savetxt('ErrorLogxArray.txt',xArray)
    np.savetxt('ErrorLogyArray.txt',yArray)
    
    
    
    
    xStitched,yStitched = Stitcher(x0,xEllipse,yEllipse,xTopStart2,yFitbutIDK,newxArr,newGeoMeanNewyArr,xArray,yArray)
    xStitchedFixed,yStitchedFixed = F.ArrayFix(xStitched,yStitched)
   
    yInterp = np.interp(xInterp,xStitchedFixed,yStitchedFixed)
   
    return yInterp

def Stitcher(x0,xElip,yElip,xCurveFit,yCurveFit,xLineFit,yLineFit,xOriginal,yOriginal):
        
    Cut = np.where(xLineFit == x0)[0][0]
        
    xLine1 = xLineFit[:Cut]
    xLine2 = xLineFit[Cut:]
    
    yLine1 = yLineFit[:Cut]
    yLine2 = yLineFit[Cut:]    
    
    LenVal1 = len(xElip)
    LenVal2 = len(xCurveFit)
    
    ActualVal1 = len(xLine1)
    ActualVal2 = len(xLine2)
    
    PadVal1 = abs(LenVal1 - ActualVal1)
    PadVal2 = abs(LenVal2 - ActualVal2)
    
    if LenVal1 < ActualVal1:
    
        xPaddedElipAlt = np.append(np.zeros(PadVal1),xElip)
        yPaddedElipAlt = np.append(np.zeros(PadVal1),yElip)
    
        MinValElipIndAlt = np.argmin(abs(yLine1 - yPaddedElipAlt))
    
        NewXElipAlt = xLine1[MinValElipIndAlt:]
        NewYElipAlt = yLine1[MinValElipIndAlt:]
        
        NewXElip2Alt = xPaddedElipAlt[:MinValElipIndAlt]
        NewYElip2Alt = yPaddedElipAlt[:MinValElipIndAlt]
    
        NewXElip2Alt = np.trim_zeros(NewXElip2Alt)
        NewYElip2Alt = np.trim_zeros(NewYElip2Alt) 
        
        FinalXInterA = np.append(NewXElip2Alt,NewXElipAlt)
        FinalYInterA = np.append(NewYElip2Alt,NewYElipAlt)
    
    else:
        xPaddedElip = np.append(np.zeros(PadVal1),xLine1)
        yPaddedElip = np.append(np.zeros(PadVal1),yLine1)
    
        MinValElipInd = np.argmin(abs(yElip-yPaddedElip))
    
        NewXElip = xPaddedElip[MinValElipInd:]
        NewYElip = yPaddedElip[MinValElipInd:]
        
        NewXElip2 = xElip[:MinValElipInd]
        NewYElip2 = yElip[:MinValElipInd]    
    
        FinalXInterA = np.append(NewXElip2,NewXElip)
        FinalYInterA = np.append(NewYElip2,NewYElip)
    
    
    
    if LenVal2 < ActualVal2:
    
        xPaddedCurveAlt = np.append(xCurveFit,np.zeros(PadVal2))
        yPaddedCurveAlt = np.append(yCurveFit,np.zeros(PadVal2))
            
        
        MinValCurveInd = np.argmin(abs(yLine2-yPaddedCurveAlt))
        
        NewXCurveAlt = xLine2[:MinValCurveInd]
        NewYCurveAlt = yLine2[:MinValCurveInd]  
        
        NewXCurve2Alt = xPaddedCurveAlt[MinValCurveInd:]
        NewYCurve2Alt = yPaddedCurveAlt[MinValCurveInd:]      
        
        NewXCurve2Alt = np.trim_zeros(NewXCurve2Alt)
        NewYCurve2Alt = np.trim_zeros(NewYCurve2Alt)   
        
        FinalXInterB = np.append(FinalXInterA,NewXCurveAlt)
        FinalXInterC = np.append(FinalXInterB,NewXCurve2Alt)
    
        FinalYInterB = np.append(FinalYInterA,NewYCurveAlt)
        FinalYInterC = np.append(FinalYInterB,NewYCurve2Alt)
        
    else: 
    
        xPaddedCurve = np.append(xLine2,np.zeros(PadVal2))
        yPaddedCurve = np.append(yLine2,np.zeros(PadVal2))
            
        
        MinValCurveInd = np.argmin(abs(yCurveFit-yPaddedCurve))
        
        NewXCurve = xPaddedCurve[:MinValCurveInd]
        NewYCurve = yPaddedCurve[:MinValCurveInd]  
        
        NewXCurve2 = xCurveFit[MinValCurveInd:]
        NewYCurve2 = yCurveFit[MinValCurveInd:]   
        
        FinalXInterB = np.append(FinalXInterA,NewXCurve)
        FinalXInterC = np.append(FinalXInterB,NewXCurve2)
    
        FinalYInterB = np.append(FinalYInterA,NewYCurve)
        FinalYInterC = np.append(FinalYInterB,NewYCurve2)
    
    
    startInds = np.where(xOriginal < min(FinalXInterC))[0]
    FinalX = np.append(xOriginal[startInds],FinalXInterC)
    FinalY = np.append(yOriginal[startInds],FinalYInterC)


    return FinalX, FinalY


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


#Disable multiplier before using
EnableCorrection = True
EnableCustomHolderCorrection = False

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
    
    
    vars()['MaxW'+str(i)] = np.array([])
    vars()['MinW'+str(i)] = np.array([])
    

    if (max(vars()['xNewMax'+str(i)]) > max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) < min(vars()['xNewMin'+str(i)])):
        for k in range(len(vars()['xNewMin'+str(i)])):
            if k == 0:
                W = 2 * (min(vars()['xNewMin'+str(i)]) - min(vars()['xNewMax'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)
    
        W = 2 * (max(vars()['xNewMax'+str(i)]) - max(vars()['xNewMin'+str(i)]))   
        vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)
        
        for j in range(1,len(vars()['xNewMax'+str(i)])):
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w) 
                
    elif (max(vars()['xNewMax'+str(i)]) < max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) < min(vars()['xNewMin'+str(i)])):
        for k in range(len(vars()['xNewMin'+str(i)])):
            if k == 0:
                W = 2 * (min(vars()['xNewMin'+str(i)]) - min(vars()['xNewMax'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)      
                
        for j in range(1,len(vars()['xNewMax'+str(i)])):
            if j == (len(vars()['xNewMax'+str(i)]) - 1):        
                w = 2 * (max(vars()['xNewMin'+str(i)]) - max(vars()['xNewMax'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
                
    elif (max(vars()['xNewMax'+str(i)]) > max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) > min(vars()['xNewMin'+str(i)])):
        for k in range(1,len(vars()['xNewMin'+str(i)])):
            if k == (len(vars()['xNewMin'+str(i)]) - 1):
                W = 2 * (max(vars()['xNewMax'+str(i)]) - max(vars()['xNewMin'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)              
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
        
        for j in range(1,len(vars()['xNewMax'+str(i)])):
            if j == (len(vars()['xNewMax'+str(i)]) - 1):
                w = 2 * (min(vars()['xNewMax'+str(i)]) - min(vars()['xNewMin'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)   
            else:
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
    
    elif (max(vars()['xNewMax'+str(i)]) < max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) > min(vars()['xNewMin'+str(i)])):
        for k in range(1,len(vars()['xNewMin'+str(i)])):
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
        
        for j in range(len(vars()['xNewMax'+str(i)])):
            if j == 0:
                w = 2 * (min(vars()['xNewMax'+str(i)]) - min(vars()['xNewMin'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)   
            else:
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w) 
            
        w = 2 * (max(vars()['xNewMin'+str(i)]) - max(vars()['xNewMax'+str(i)]))   
        vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
        
        
    vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)] + ((Slitwidth*vars()['yNewMaxUnCorr'+str(i)])/vars()['MaxW'+str(i)])**2
    vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)] + ((Slitwidth*vars()['yNewMinUnCorr'+str(i)])/vars()['MinW'+str(i)])**2
 

    # vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)]
    # vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)]

    vars()['yPMaxAlt'+str(i)] = Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)],0.6)
    vars()['yPMinAlt'+str(i)] = Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)],0.6)


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
    plt.plot(xP,vars()['yPMaxAlt'+str(i)], color = 'grey', linestyle="dotted", label = "TM2")
    plt.plot(xP,vars()['yPMinAlt'+str(i)], color = 'orange', linestyle="dotted", label = "Tm2")    
    
    plt.plot(xP,vars()['yP'+files1[0]], color = 'orange', label = "Substrate")
    
    plt.legend()
    
