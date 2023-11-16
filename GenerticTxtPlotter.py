import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import Functions1 as F1
import nPlotterFunctions as nPF
import csv
from scipy.optimize import curve_fit


path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\FitOptTest\\raw'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))

path2, dirs2, files2, data2, header2 = nPF.fileGrabber('GenerticTxtPlotter.py','\\FitOptTest\\csv','.csv')

#In nm
Slitwidth = 2

xMin = 200
xMax = 800

sizeFigSqaure = F1.CTI(20)

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
    
    #vars()[files1[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files1[i]],xMin,xMax))
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

for i in range(len(files2)):
    vars()[f'xVals {files2[i]}'] = np.array(F1.ListExtract(data2[i],0))
    vars()[f'yVals {files2[i]}'] = np.array(F1.ListExtract(data2[i],1))

#Since we know the cutoff is 0.6
#Given we know how the array is constructed we can assume the next index is for the first maxima/minima
for i in range(len(files2)):
    vars()[f'yInd {files2[i]}'] = np.where(vars()[f'yVals {files2[i]}'] == F1.FindNearestVal(vars()[f'yVals {files2[i]}'],0.6))[0][0]
    vars()[f'y1 {files2[i]}'] = vars()[f'yVals {files2[i]}'][vars()[f'yInd {files2[i]}'] - 4]
    vars()[f'x1 {files2[i]}'] = vars()[f'xVals {files2[i]}'][vars()[f'yInd {files2[i]}'] - 4]
    vars()[f'y2 {files2[i]}'] = vars()[f'yVals {files2[i]}'][vars()[f'yInd {files2[i]}'] + 1]
    vars()[f'x2 {files2[i]}'] = vars()[f'xVals {files2[i]}'][vars()[f'yInd {files2[i]}'] + 1]
    vars()[f'x0 {files2[i]}'], vars()[f'y0 {files2[i]}'], vars()[f'a {files2[i]}'], vars()[f'b {files2[i]}'], vars()[f'c {files2[i]}'] = F1.ellipseParamFinder(vars()[f'x1 {files2[i]}'],vars()[f'y1 {files2[i]}'],vars()[f'x2 {files2[i]}'],vars()[f'y2 {files2[i]}'])     
    vars()[f'xEllipse {files2[i]}'] = np.arange(vars()[f'x1 {files2[i]}'], vars()[f'x2 {files2[i]}']+1, 1)
    vars()[f'yEllipse {files2[i]}'] = F1.quarterEllipseFinder(vars()[f'xEllipse {files2[i]}'],vars()[f'a {files2[i]}'],vars()[f'b {files2[i]}'], vars()[f'x0 {files2[i]}'],vars()[f'y0 {files2[i]}'])

guess = np.array([np.exp(1),1,260])

for i in range(len(files2)):
    vars()[f'yTopStart {files2[i]}'] =  vars()[f'yVals {files2[i]}'][vars()[f'yInd {files2[i]}']+1:]
    vars()[f'xTopStart {files2[i]}'] =  vars()[f'xVals {files2[i]}'][vars()[f'yInd {files2[i]}']+1:]
    
    
    vars()[f'TransFit {files2[i]}'], vars()[f'Resi {files2[i]}'] = curve_fit(F1.logCurve,vars()[f'xTopStart {files2[i]}'],vars()[f'yTopStart {files2[i]}'],guess, maxfev=500000000)
       
    vars()[f'xTopStart2 {files2[i]}'] = np.arange(vars()[f'xVals {files2[i]}'][vars()[f'yInd {files2[i]}']+1], 801, 1)
        
    vars()[f'yFitbutIDK {files2[i]}'] = F1.logCurve(vars()[f'xTopStart2 {files2[i]}'],vars()[f'TransFit {files2[i]}'][0],vars()[f'TransFit {files2[i]}'][1],vars()[f'TransFit {files2[i]}'][2])





#Since step size is 1 we can use np.gradient

#!!! Note: Change X values to cover all x points for final fit

for i in range(len(files2)):
    vars()[f'gradEllipse {files2[i]}'] = np.gradient(vars()[f'yEllipse {files2[i]}'])
    vars()[f'bEllipse {files2[i]}'] = vars()[f'yEllipse {files2[i]}'] - (vars()[f'gradEllipse {files2[i]}']*vars()[f'xEllipse {files2[i]}'])
    
    vars()[f'gradyFitbutIDK {files2[i]}'] = np.gradient(vars()[f'yFitbutIDK {files2[i]}'])
    vars()[f'bFitbutIDK {files2[i]}'] = vars()[f'yFitbutIDK {files2[i]}'] - (vars()[f'gradyFitbutIDK {files2[i]}']*vars()[f'xTopStart2 {files2[i]}'])
 
allowance = 20
for i in range(len(files2)):
    vars()[f'xArr {files2[i]}'] = np.linspace(vars()[f'x0 {files2[i]}']-allowance,vars()[f'x0 {files2[i]}']+allowance,1001)
    vars()[f'IncLoc {files2[i]}'] = F1.LineComp(vars()[f'xArr {files2[i]}'],vars()[f'gradEllipse {files2[i]}'],vars()[f'gradyFitbutIDK {files2[i]}'],vars()[f'bEllipse {files2[i]}'],vars()[f'bFitbutIDK {files2[i]}'],5)


for i in range(len(files2)):
    vars()[f'yArr {files2[i]}'] = (vars()[f'xArr {files2[i]}']*vars()[f'gradEllipse {files2[i]}'][vars()[f'IncLoc {files2[i]}'][1]] + vars()[f'bEllipse {files2[i]}'][vars()[f'IncLoc {files2[i]}'][1]])
    vars()[f'yArr1 {files2[i]}'] = (vars()[f'xArr {files2[i]}']*vars()[f'gradyFitbutIDK {files2[i]}'][vars()[f'IncLoc {files2[i]}'][2]] + vars()[f'bFitbutIDK {files2[i]}'][vars()[f'IncLoc {files2[i]}'][2]])
    
    vars()[f'xArr {files2[i]}'] = np.round(vars()[f'xArr {files2[i]}']).astype(int)

    vars()[f'GeoMeanNewyArr {files2[i]}'] = np.sqrt((vars()[f'yArr {files2[i]}'] * vars()[f'yArr1 {files2[i]}']))
    
for i in range(len(files2)):
    vars()[f'newxArr {files2[i]}'], vars()[f'newGeoMeanNewyArr {files2[i]}'] = F1.ArrayCondenser(vars()[f'xArr {files2[i]}'],vars()[f'GeoMeanNewyArr {files2[i]}'])


for i in range(len(files2)):
    vars()[f'xStitched {files2[i]}'],vars()[f'yStitched {files2[i]}'] = F1.Stitcher(vars()[f'x0 {files2[i]}'],vars()[f'xEllipse {files2[i]}'],vars()[f'yEllipse {files2[i]}'],vars()[f'xTopStart2 {files2[i]}'],vars()[f'yFitbutIDK {files2[i]}'],vars()[f'newxArr {files2[i]}'],vars()[f'newGeoMeanNewyArr {files2[i]}'],vars()[files[0]+'T'][0],vars()['yFiltered'+str(0)])
    vars()[f'xStitchedFixed {files2[i]}'],vars()[f'yStitchedFixed {files2[i]}'] = F1.ArrayFix(vars()[f'xStitched {files2[i]}'],vars()[f'yStitched {files2[i]}'])
boolarray = [True,False]

for i in range(len(files2)):
    vars()[f'yStichedScorrected {files2[i]}'] = F1.fitCorrectingFunction(vars()['yFiltered'+str(0)],vars()[f'yStitchedFixed {files2[i]}'],boolarray[i])

# plt.figure(0, figsize = (11.7,8.3), dpi = 300)
plt.figure(0)
# plt.xlim(225,400)
plt.xlim(250,500)
plt.ylim(0.6,0.85)
for i in range(len(files2)):
    plt.scatter(vars()[f'xVals {files2[i]}'],vars()[f'yVals {files2[i]}'])
    plt.scatter(vars()[f'x0 {files2[i]}'], vars()[f'y0 {files2[i]}'])
    # plt.plot(vars()[f'xTopStart2 {files2[i]}'],vars()[f'yFitbutIDK {files2[i]}'])
    # plt.plot(vars()[f'xEllipse {files2[i]}'],vars()[f'yEllipse {files2[i]}'])
    # plt.plot(vars()[f'newxArr {files2[i]}'],vars()[f'newGeoMeanNewyArr {files2[i]}'])
    plt.plot(vars()[f'xStitchedFixed {files2[i]}'],vars()[f'yStichedScorrected {files2[i]}'])
plt.plot(vars()[files[0]+'T'][0],vars()['yFiltered'+str(0)])


#(np.round(vars()[f'xArr {files2[i]}'])).astype(int)

# for i in range(len(files)):            
#     vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])
    
    
#     #0.075 is standard value for the second variable in the function. but adjust as neccessary
#     vars()['xNewMax'+str(i)], vars()['yNewMaxUnCorr'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMinUnCorr'+str(i)] = F.BackNodeFix(3,0.075,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)])
    
#     vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)]
#     vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)]




# # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# plt.figure(0,figsize = (8,5),dpi = 600)
# plt.minorticks_on()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# for i in range(len(files)):
#     plt.plot(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)], label = legend[i])
# plt.title('Comparison of mounting techniques')
# plt.xlabel("Wavelength (nm)")
# plt.ylabel("Transmission/Absorbance")
 
# # plt.plot(xP,vars()['yP'+files1[0]], color = 'orange', label = "Substrate")

# plt.legend()

