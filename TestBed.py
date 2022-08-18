import numpy as np
import matplotlib.pyplot as plt

def gauss(mean,std,array):
   yarray = np.zeros(len(array)) 
   for x in range(0,len(array)):
        yarray[x] = (1 / (std * np.sqrt(2*np.pi))) * np.exp(-((array[x]-mean)**2)/(2*std**2))
   return yarray

Lambda = np.array([292,302,314,327,341,358,376,397,421,449,482,519,565,618,685,769,296,307,319,333,349,366,386,408,434,464,499,540,590,649,722])
Thicc = np.array([4635.34456,2165.840875,2683.106913,1419.97695,1611.031029,1746.734817,1677.946119,1591.554847,1607.481509,1595.678189,1561.901806,1364.72997,1791.615589,1590.266815,1776.040617,1872.726431,2658.374403,4116.196434,1236.028022,1951.660233,1571.435833,1776.551548,1619.634968,1590.113786,1686.409821,1459.758653,1612.8373,1339.562859,1920.921311,1583.824594,1833.15318])

ThiccInds = Thicc.argsort()
LambdaSort = Lambda[ThiccInds[::-1]]
ThiccSort = Thicc[ThiccInds[::-1]]


ControlCutOff = 2

ThiccSort = ThiccSort[ControlCutOff:]
LambdaSort = LambdaSort[ControlCutOff:]

# TargetVal = (np.mean(ThiccSort)) / 10

Std = np.std(ThiccSort)
Mean = np.mean(ThiccSort)

ThiccGauss = gauss(Mean,Std,ThiccSort)

Multi = max(ThiccGauss) / max(LambdaSort)

plt.figure(1)
plt.scatter(ThiccSort,ThiccGauss)

plt.scatter(ThiccSort,Multi * LambdaSort)

def ThicknessScrub(Lambda,d,CutOff):
    
    dInds = d.argsort()
    d = d[dInds[::-1]]
    Lambda = Lambda[dInds[::-1]]
    
    d = d[CutOff:]
    Lambda = Lambda[CutOff:]
    
    dStd = np.std(d)
    dMean = np.mean(d)
    
    dGauss = gauss(dMean,dStd,d)
    
    Multiplier = max(dGauss) / max(Lambda)
    
    GoodInds = np.where(Multiplier*Lambda < dGauss)
    BadInds = np.where(Multi*LambdaSort >= ThiccGauss)
    
                
    RejectedLambda = Lambda[BadInds]
    Newd = d[GoodInds]
    newMean = np.mean(Newd)
    error = np.std(Newd) / np.sqrt(len(Newd))
    
    
    return newMean,error,RejectedLambda