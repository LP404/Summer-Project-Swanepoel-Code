import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import nPlotterFunctions as nPF
import csv
from scipy.optimize import curve_fit
import copy

def FindNearestVal(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

def CTI(value):
    return value/2.54

def ArrayCondenser(x,y):
      newXarray = np.array([])
      newYarray = np.array([])
    
      newXarray = np.append(newXarray,x[0])
      backTotalCount = 0
      backTotal = 0
    
      for i in range(1,len(x)):
          if x[i] == x[i-1]:
             backTotalCount += 1
          else:   
             newXarray = np.append(newXarray,x[i])
             newYarray = np.append(newYarray,np.mean(y[backTotal:i]))
             backTotal = backTotalCount + 1
    
    
      backTotal = backTotalCount + 1
      newYarray = np.append(newYarray,np.mean(y[backTotal:i]))
    
      return newXarray,newYarray

def ArrayFix(x,y):
    
    x2 = np.array([x[0]])
    y2 = np.array([y[0]])
    
    for i in range(1,len(x)):
        if x[i] != x[i-1]:
            x2 = np.append(x2,x[i])
            y2 = np.append(y2,y[i])
        else:
            yHold = (y[i] + y[i-1]) / 2
            y2[i-1] = yHold
    
    return x2, y2

def LineIntersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

#For a horizonal ellipse
def ellipseParamFinder(x1,y1,x2,y2):
    line1 = [[200,y1] , [800,y1]]
    line2 = [[x2,0] , [x2,1]]
    x0, y0 = LineIntersection(line1 ,line2)
    a = abs(x1-x0)
    b = abs(y2-y0)
    c = np.sqrt((a**2 + b**2))
    
    return x0,y0,a,b,c

def quarterEllipseFinder(x,a,b,x0,y0):
    y = y0 + (((b/a) * (np.sqrt((a**2-((x-x0)**2))))))
    return y
    
def logCurve(x,a,b,c):
    
    y = np.emath.logn(a,x-c) + b
        
    return y 

#Note define an x before passing in
def LineComp(x,m1,m2,c1,c2,delta):
    m1Index = 0
    m2Index = 0
    MaxOccur = 0
    for i in range(len(m1)):
        y1 = m1[i]*x + c1[i]
        for j in range(len(m2)):
            y2 = m2[j]*x + c2[j]
            checkArray = y1/y2
            for k in range(len(checkArray)):
                if checkArray[k] < 1:
                    checkArray[k] = (1 / checkArray[k])
                else:
                    pass
                
            oneOccur = np.count_nonzero(np.around(checkArray,delta) == 1)
            if oneOccur > MaxOccur:
                m1Index = i
                m2Index = j
                MaxOccur = oneOccur
            else:
                pass
    return MaxOccur, m1Index, m2Index

    
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
    
    xPaddedElip = np.append(np.zeros(PadVal1),xLine1)
    yPaddedElip = np.append(np.zeros(PadVal1),yLine1)
    
    xPaddedCurve = np.append(xLine2,np.zeros(PadVal2))
    yPaddedCurve = np.append(yLine2,np.zeros(PadVal2))
    
    MinValElipInd = np.argmin(abs(yElip-yPaddedElip))
    MinValCurveInd = np.argmin(abs(yCurveFit-yPaddedCurve))
    
    NewXElip = xPaddedElip[MinValElipInd:]
    NewYElip = yPaddedElip[MinValElipInd:]
    
    NewXElip2 = xElip[:MinValElipInd]
    NewYElip2 = yElip[:MinValElipInd]
    
    NewXCurve = xPaddedCurve[:MinValCurveInd]
    NewYCurve = yPaddedCurve[:MinValCurveInd]  
    
    NewXCurve2 = xCurveFit[MinValCurveInd:]
    NewYCurve2 = yCurveFit[MinValCurveInd:]    
    
    FinalXInterA = np.append(NewXElip2,NewXElip)
    FinalXInterB = np.append(FinalXInterA,NewXCurve)
    FinalXInterC = np.append(FinalXInterB,NewXCurve2)
    
    startInds = np.where(xOriginal < min(FinalXInterC))[0]
    FinalX = np.append(xOriginal[startInds],FinalXInterC)
    
    FinalYInterA = np.append(NewYElip2,NewYElip)
    FinalYInterB = np.append(FinalYInterA,NewYCurve)
    FinalYInterC = np.append(FinalYInterB,NewYCurve2)
    
    FinalY = np.append(yOriginal[startInds],FinalYInterC)
        
    return FinalX, FinalY
    
def fitCorrectingFunction(y,yFit,MaxMin):
    
    yFitDup = copy.deepcopy(yFit)

    

    for i in range(len(y)):
        
        if yFitDup[i] < y[i] and MaxMin == True:
        
            corrFactor = y[i]/yFitDup[i]
                    
            yFitDup[i] = corrFactor*yFitDup[i]
            
        
        elif yFitDup[i] > y[i] and MaxMin == False:
        
            corrFactor = y[i]/yFitDup[i]
            
            yFitDup[i] = corrFactor*yFitDup[i]
            
            
    return yFitDup

def main():
    return

if __name__ == '__main__':
    main()
