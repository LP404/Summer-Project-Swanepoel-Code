import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import Functions1 as F1
import nPlotterFunctions as nPF
import csv
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi
from natsort import natsorted
import copy


def CTI(value):
    return value/2.54

def Normalise(x):
    
    xMax = max(x)
    xMin = min(x)
    xDiff = xMax - xMin
    
    xNorm = np.zeros_like(x)
    
    for i in range(len(x)):
        xNorm[i] = (x[i] - xMin) / xDiff
    
    return xNorm

def Trim(Array,start,end):
    delPoints = np.where((Array[0] < start) | (Array[0] > end))[0]
    
    newArray = np.array([np.delete(Array[0],delPoints),np.delete(Array[1],delPoints)])
    
    return newArray


def PlankLaw(T,x):
    
    h = 6.62607015e-34
    c = 299792458
    k = 1.380649e-23
    x = x*1e-9
    x = c/x
    
    A1 = (2*h*(x**3)) / (c**2)
    expo = ((h*x) / (k*T))
    A2 = 1 / (np.exp(expo) - 1)
    
    return A1*A2


def ConvertDependance(x,y):
    
    h = 6.62607015e-34
    c = 299792458
    x = x*1e-9
    
    xNew = ((h*c) / (x)) * 6.242e18
    yNew = y * (((x)**2)/(h*c))

    return xNew,yNew

def NoiseFilter(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    yFiltered = filtfilt(b, a, yArray)
    
    return yFiltered

xMin = 280
xMax = 700

sizeFigSqaure = F1.CTI(20)

pathA, dirsA, fileNamesA = next(os.walk(os.path.dirname(os.path.realpath('GaussFitPlotter.py')) + '\\fitykExports\\TempGaussFits\\DataExportAlpha\\NewQ'))
pathB, dirsB, fileNamesB = next(os.walk(os.path.dirname(os.path.realpath('GaussFitPlotter.py')) + '\\fitykExports\\TempGaussFits\\DataExportBeta\\NewQ'))
pathK, dirsK, fileNamesK = next(os.walk(os.path.dirname(os.path.realpath('GaussFitPlotter.py')) + '\\fitykExports\\TempGaussFits\\DataExportKappa\\NewQ'))

# pathTu, dirsTu, fileNamesTu = next(os.walk(os.path.dirname(os.path.realpath('GaussFitPlotter.py')) + '\\fitykExports\\TempGaussFits\\Tumadir'))


fileNamesA = natsorted(fileNamesA)
fileNamesB = natsorted(fileNamesB)
fileNamesK = natsorted(fileNamesK)
# fileNamesTu = natsorted(fileNamesTu)

TempValsA = np.array([])
TempValsB = np.array([])
TempValsK = np.array([])
# TempValsTu = np.array([])


for i in range(len(fileNamesA)):
    fileNamesA[i] = fileNamesA[i][:-4]
    vars()[fileNamesA[i] + ' header alpha'], vars()[fileNamesA[i] + ' data alpha'] = nPF.datGrabberAlt(fileNamesA[i],"fitykExports\\TempGaussFits\\DataExportAlpha\\NewQ\\")
    
    PeakArrays = copy.deepcopy(vars()[fileNamesA[i] + ' data alpha'][3: 3 + vars()[fileNamesA[i] + ' header alpha'][-1]])
    
    InitList = []
    
    for j in range(len(PeakArrays)):
        
        Loc = np.where(vars()[fileNamesA[i] + ' data alpha'][3+j] == max(PeakArrays[j]))[0][0]
        
        InitList.append(vars()[fileNamesA[i] + ' data alpha'][0][Loc])
        
    InitCheck = [[InitList[k],k,PeakArrays[k]] for k in range(len(InitList))]
    InitCheck.sort(key=lambda InitCheck: InitCheck[0])
    
    vars()[fileNamesA[i] + ' IntInt alpha'] = np.array([])
    
    for j in range(len(PeakArrays)):
        vars()[fileNamesA[i] + ' data alpha'][3 + j] = InitCheck[j][2]
        
        RunTot = 0
        for k in range(len(InitCheck[j][2])):
            RunTot += (InitCheck[j][2][k] - max(InitCheck[j][2]))**2
        Var = RunTot/len(InitCheck[j][2])
        Std = Var**0.5
        IntInt = max(InitCheck[j][2]) * 2*np.pi * Std
        
        vars()[fileNamesA[i] + ' IntInt alpha'] = np.append(vars()[fileNamesA[i] + ' IntInt alpha'],IntInt)
        
        
    TempValsA = np.append(TempValsA,float(fileNamesA[i][:-1]))
    
    vars()[fileNamesA[i] + ' peak locs alpha'] = np.array([])
    
    for j in range(len(vars()[fileNamesA[i] + ' header alpha'][3:3 + vars()[fileNamesA[i] + ' header alpha'][-1]])):
        PeakLoc = vars()[fileNamesA[i] + ' data alpha'][0][np.where(max(vars()[fileNamesA[i] + ' data alpha'][3+j]) == (vars()[fileNamesA[i] + ' data alpha'][3+j]))[0][0]]
        
        vars()[fileNamesA[i] + ' peak locs alpha'] = np.append(vars()[fileNamesA[i] + ' peak locs alpha'],PeakLoc)


for i in range(len(fileNamesB)):
    fileNamesB[i] = fileNamesB[i][:-4]
    vars()[fileNamesB[i] + ' header beta'], vars()[fileNamesB[i] + ' data beta'] = nPF.datGrabberAlt(fileNamesB[i],"fitykExports\\TempGaussFits\\DataExportBeta\\NewQ\\")

    PeakArrays = copy.deepcopy(vars()[fileNamesB[i] + ' data beta'][3: 3 + vars()[fileNamesB[i] + ' header beta'][-1]])
    
    InitList = []
    
    for j in range(len(PeakArrays)):
        
        Loc = np.where(vars()[fileNamesB[i] + ' data beta'][3+j] == max(PeakArrays[j]))[0][0]
        
        InitList.append(vars()[fileNamesB[i] + ' data beta'][0][Loc])
        
    InitCheck = [[InitList[k],k,PeakArrays[k]] for k in range(len(InitList))]
    InitCheck.sort(key=lambda InitCheck: InitCheck[0])
    
    vars()[fileNamesA[i] + ' IntInt beta'] = np.array([])
    
    for j in range(len(PeakArrays)):
        vars()[fileNamesB[i] + ' data beta'][3 + j] = InitCheck[j][2]

        RunTot = 0
        for k in range(len(InitCheck[j][2])):
            RunTot += (InitCheck[j][2][k] - max(InitCheck[j][2]))**2
        Var = RunTot/len(InitCheck[j][2])
        Std = Var**0.5
        IntInt = max(InitCheck[j][2]) * 2*np.pi * Std
        
        vars()[fileNamesA[i] + ' IntInt beta'] = np.append(vars()[fileNamesA[i] + ' IntInt beta'],IntInt)


    TempValsB = np.append(TempValsB,float(fileNamesB[i][:-1]))

    vars()[fileNamesB[i] + ' peak locs beta'] = np.array([])
    
    for j in range(len(vars()[fileNamesB[i] + ' header beta'][3:3 + vars()[fileNamesB[i] + ' header beta'][-1]])):
        PeakLoc = vars()[fileNamesB[i] + ' data beta'][0][np.where(max(vars()[fileNamesB[i] + ' data beta'][3+j]) == (vars()[fileNamesB[i] + ' data beta'][3+j]))[0][0]]
        
        vars()[fileNamesB[i] + ' peak locs beta'] = np.append(vars()[fileNamesB[i] + ' peak locs beta'],PeakLoc)


for i in range(len(fileNamesK)):
    fileNamesK[i] = fileNamesK[i][:-4]
    vars()[fileNamesK[i] + ' header kappa'], vars()[fileNamesK[i] + ' data kappa'] = nPF.datGrabberAlt(fileNamesK[i],"fitykExports\\TempGaussFits\\DataExportKappa\\NewQ\\")

    PeakArrays = copy.deepcopy(vars()[fileNamesK[i] + ' data kappa'][3: 3 + vars()[fileNamesK[i] + ' header kappa'][-1]])
    
    InitList = []
    
    for j in range(len(PeakArrays)):
        
        Loc = np.where(vars()[fileNamesK[i] + ' data kappa'][3+j] == max(PeakArrays[j]))[0][0]
        
        InitList.append(vars()[fileNamesK[i] + ' data kappa'][0][Loc])
        
    InitCheck = [[InitList[k],k,PeakArrays[k]] for k in range(len(InitList))]
    InitCheck.sort(key=lambda InitCheck: InitCheck[0])
    
    vars()[fileNamesA[i] + ' IntInt kappa'] = np.array([])
    
    for j in range(len(PeakArrays)):
        vars()[fileNamesK[i] + ' data kappa'][3 + j] = InitCheck[j][2]

        RunTot = 0
        for k in range(len(InitCheck[j][2])):
            RunTot += (InitCheck[j][2][k] - max(InitCheck[j][2]))**2
        Var = RunTot/len(InitCheck[j][2])
        Std = Var**0.5
        IntInt = max(InitCheck[j][2]) * 2*np.pi * Std
        
        vars()[fileNamesA[i] + ' IntInt kappa'] = np.append(vars()[fileNamesA[i] + ' IntInt kappa'],IntInt)


    TempValsK = np.append(TempValsK,float(fileNamesK[i][:-1]))
    
    vars()[fileNamesK[i] + ' peak locs kappa'] = np.array([])
    
    for j in range(len(vars()[fileNamesK[i] + ' header kappa'][3:3 + vars()[fileNamesK[i] + ' header kappa'][-1]])):
        PeakLoc = vars()[fileNamesK[i] + ' data kappa'][0][np.where(max(vars()[fileNamesK[i] + ' data kappa'][3+j]) == (vars()[fileNamesK[i] + ' data kappa'][3+j]))[0][0]]
        
        vars()[fileNamesK[i] + ' peak locs kappa'] = np.append(vars()[fileNamesK[i] + ' peak locs kappa'],PeakLoc)




# for i in range(len(fileNamesTu)):
#     fileNamesTu[i] = fileNamesTu[i][:-4]
#     vars()[fileNamesTu[i] + ' header Tu'], vars()[fileNamesTu[i] + ' data Tu'] = nPF.datGrabberAlt(fileNamesTu[i],"fitykExports\\TempGaussFits\\Tumadir\\")

#     TempValsTu = np.append(TempValsTu,float(fileNamesTu[i][11:-9])) 
    
#     vars()[fileNamesTu[i] + ' peak locs Tu'] = np.array([])
    
#     for j in range(len(vars()[fileNamesTu[i] + ' header Tu'][3:3 + vars()[fileNamesTu[i] + ' header Tu'][-1]])):
#         PeakLoc = vars()[fileNamesTu[i] + ' data Tu'][0][np.where(max(vars()[fileNamesTu[i] + ' data Tu'][3+j]) == (vars()[fileNamesTu[i] + ' data Tu'][3+j]))[0][0]]
        
#         vars()[fileNamesTu[i] + ' peak locs Tu'] = np.append(vars()[fileNamesTu[i] + ' peak locs Tu'],PeakLoc)
        
    # vars()[fileNamesTu[i] + ' peak locs Tu'].sort()

headerAlphaPL, dataAlphaPL = nPF.datGrabberAlt('300K','fitykExports\\PLdata\\a\\')
headerAlphaPLlowT, dataAlphaPLlowT = nPF.datGrabberAlt('23K','fitykExports\\PLdata\\a\\')

headerBetaPL, dataBetaPL = nPF.datGrabberAlt('300K','fitykExports\\PLdata\\b\\')
headerBetaPLlowT, dataBetaPLlowT = nPF.datGrabberAlt('23K','fitykExports\\PLdata\\b\\')

headerKappaPL, dataKappaPL = nPF.datGrabberAlt('300K','fitykExports\\PLdata\\k\\')
headerKappaPLlowT, dataKappaPLlowT = nPF.datGrabberAlt('23K','fitykExports\\PLdata\\k\\')

GaussColours = ['red','green','dodgerblue','darkblue','rebeccapurple','darkorchid','darkmagenta']


for i in range(len(TempValsA)):
    plt.figure(0, dpi = 300)
    plt.xlabel('1/T (K)')
    plt.ylabel('Peak energy (eV)')
    plt.title('Peak position variation of alpha')
    plt.scatter(1/TempValsA[i],vars()[fileNamesA[i] + ' peak locs alpha'][0], color = GaussColours[2])
    plt.scatter(1/TempValsA[i],vars()[fileNamesA[i] + ' peak locs alpha'][1], color = GaussColours[3])
    plt.scatter(1/TempValsA[i],vars()[fileNamesA[i] + ' peak locs alpha'][2], color = GaussColours[5])
    plt.scatter(1/TempValsA[i],vars()[fileNamesA[i] + ' peak locs alpha'][3], color = GaussColours[6])

    plt.figure(1, dpi = 300)
    plt.xlabel('T (K)')
    plt.ylabel('Integrated Intensity (Arb)')
    plt.title('Integrated instensity variation of alpha (log)')
    plt.yscale('log')
    plt.scatter(TempValsA[i],vars()[fileNamesA[i] + ' IntInt alpha'][0], color = GaussColours[2])
    plt.scatter(TempValsA[i],vars()[fileNamesA[i] + ' IntInt alpha'][1], color = GaussColours[3])
    plt.scatter(TempValsA[i],vars()[fileNamesA[i] + ' IntInt alpha'][2], color = GaussColours[5])
    plt.scatter(TempValsA[i],vars()[fileNamesA[i] + ' IntInt alpha'][3], color = GaussColours[6])
 

for i in range(len(TempValsB)):
    plt.figure(2, dpi = 300)
    plt.xlabel('1/T (K)')
    plt.ylabel('Peak energy (eV)')
    plt.title('Peak position variation of beta')
    plt.scatter(1/TempValsB[i],vars()[fileNamesB[i] + ' peak locs beta'][0], color = GaussColours[1])
    plt.scatter(1/TempValsB[i],vars()[fileNamesB[i] + ' peak locs beta'][1], color = GaussColours[2])
    plt.scatter(1/TempValsB[i],vars()[fileNamesB[i] + ' peak locs beta'][2], color = GaussColours[5])
    plt.scatter(1/TempValsB[i],vars()[fileNamesB[i] + ' peak locs beta'][3], color = GaussColours[6])
    
    plt.figure(3, dpi = 300)
    plt.xlabel('T (K)')
    plt.ylabel('Integrated Intensity (Arb)')
    plt.title('Integrated instensity variation of beta (log)')
    plt.yscale('log')
    plt.scatter(TempValsB[i],vars()[fileNamesB[i] + ' IntInt beta'][0], color = GaussColours[1])
    plt.scatter(TempValsB[i],vars()[fileNamesB[i] + ' IntInt beta'][1], color = GaussColours[2])
    plt.scatter(TempValsB[i],vars()[fileNamesB[i] + ' IntInt beta'][2], color = GaussColours[5])
    plt.scatter(TempValsB[i],vars()[fileNamesB[i] + ' IntInt beta'][3], color = GaussColours[6])
    
for i in range(len(TempValsK)):
    plt.figure(4, dpi = 300)
    plt.xlabel('1/T (K)')
    plt.ylabel('Peak energy (eV)')
    plt.title('Peak position variation of kappa')
    plt.scatter(1/TempValsK[i],vars()[fileNamesK[i] + ' peak locs kappa'][0], color = 'orange')
    plt.scatter(1/TempValsK[i],vars()[fileNamesK[i] + ' peak locs kappa'][1], color = GaussColours[2])
    plt.scatter(1/TempValsK[i],vars()[fileNamesK[i] + ' peak locs kappa'][2], color = GaussColours[5])
    plt.scatter(1/TempValsK[i],vars()[fileNamesK[i] + ' peak locs kappa'][3], color = GaussColours[6])
    
    plt.figure(5, dpi = 300)
    plt.xlabel('T (K)')
    plt.ylabel('Integrated Intensity (Arb)')
    plt.title('Integrated instensity variation of kappa (log)')
    plt.yscale('log')
    plt.scatter(TempValsK[i],vars()[fileNamesK[i] + ' IntInt kappa'][0], color = 'orange')
    plt.scatter(TempValsK[i],vars()[fileNamesK[i] + ' IntInt kappa'][1], color = GaussColours[2])
    plt.scatter(TempValsK[i],vars()[fileNamesK[i] + ' IntInt kappa'][2], color = GaussColours[5])
    plt.scatter(TempValsK[i],vars()[fileNamesK[i] + ' IntInt kappa'][3], color = GaussColours[6])


# headerAlphaPL, dataAlphaPL = nPF.datGrabberSimple('300K','fitykExports\\PLdata\\a\\')
# headerAlphaPLlowT, dataAlphaPLlowT = nPF.datGrabberSimple('23K','fitykExports\\PLdata\\a\\')

# headerBetaPL, dataBetaPL = nPF.datGrabberSimple('300K','fitykExports\\PLdata\\b\\')
# headerBetaPLlowT, dataBetaPLlowT = nPF.datGrabberSimple('23K','fitykExports\\PLdata\\b\\')

# headerKappaPL, dataKappaPL = nPF.datGrabberSimple('300K','fitykExports\\PLdata\\k\\')
# headerKappaPLlowT, dataKappaPLlowT = nPF.datGrabberSimple('23K','fitykExports\\PLdata\\k\\')



# Red = np.linspace(0,255,len(files)+1).astype('int32')
# Blue = np.flip(Red)


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple']
LineStyle = ['solid','dashed']

TickSize = 14
LabelSize = 14
LegSize = 10


GraphResolution = 100
sizeFigSqaure = CTI(30)

# plt.figure(83, figsize = (sizeFigSqaure * 2, sizeFigSqaure),dpi = GraphResolution)
# plt.minorticks_on()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')

# plt.plot(dataAlphaCLLT[0],dataAlphaCLLT[1]/max(dataAlphaCLLT[1]), label = LabList[0] + ' 80K', color = LineColour[0],linestyle = LineStyle[0])
# plt.plot(dataBetaCLLT[0],dataBetaCLLT[1]/max(dataBetaCLLT[1]), label = LabList[1] + ' 80K', color = LineColour[1],linestyle = LineStyle[0])
# plt.plot(dataKappaCLLT[0],dataKappaCLLT[1]/max(dataKappaCLLT[1]), label = LabList[2] + ' 80K', color = LineColour[2],linestyle = LineStyle[0])

# plt.plot(dataAlphaCL[0],dataAlphaCL[1]/max(dataAlphaCL[1]), label = LabList[0] + ' 293K', color = LineColour[0],linestyle = LineStyle[1])
# plt.plot(dataBetaCL[0],dataBetaCL[1]/max(dataBetaCL[1]), label = LabList[1] + ' 293K', color = LineColour[1],linestyle = LineStyle[1])
# plt.plot(dataKappaCL[0],dataKappaCL[1]/max(dataKappaCL[1]), label = LabList[2] + ' 293K', color = LineColour[2],linestyle = LineStyle[1])

# plt.xlim(1.78622,4.5)
# plt.ylim(0,1.05)

# plt.title(f'Cahtodoluminescence of Ga₂O₃ polymorphs', fontsize = 18)
# plt.xlabel("Energy (eV)", fontsize = 18)
# plt.ylabel("Intensity (counts/eV)", fontsize = 18)
# plt.legend(loc = 'best', fontsize = 26)

GaussColours = ['red','green','dodgerblue','darkblue','rebeccapurple','darkorchid','darkmagenta']

GraphResolution = 300
sizeFigSqaure = CTI(30)

# plt.figure(0,dpi = GraphResolution)
# plt.minorticks_on()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# plt.yscale('log')
# # plt.xlim(1.78622,4.5)
# plt.ylim(2e-3,1)

# plt.xlabel("Energy (eV)", fontsize = 12)
# plt.ylabel("Intensity (Arb. units)", fontsize = 12)

# plt.scatter(dataAlphaPL[0],dataAlphaPL[1]/max(dataAlphaPL[1]), label = LabList[0] + ' 23K', color = 'darkcyan', s = 1)
# plt.plot(dataAlphaPL[0],dataAlphaPL[2]/max(dataAlphaPL[1]), label = '3.9eV', color = GaussColours[6])
# plt.plot(dataAlphaPL[0],dataAlphaPL[3]/max(dataAlphaPL[1]), label = '3.5eV', color = GaussColours[5])
# plt.plot(dataAlphaPL[0],dataAlphaPL[4]/max(dataAlphaPL[1]), label = '3.1eV', color = GaussColours[3])
# plt.plot(dataAlphaPL[0],dataAlphaPL[5]/max(dataAlphaPL[1]), label = '2.8eV', color = GaussColours[2])
# plt.plot(dataAlphaPL[0],dataAlphaPL[6]/max(dataAlphaPL[1]), color = 'k', linestyle = '--')
# plt.legend(loc = 'best', fontsize = 12)

# plt.figure(1,dpi = GraphResolution)
# plt.minorticks_on()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# plt.yscale('log')
# # plt.xlim(1.78622,4.5)
# plt.ylim(2e-3,1)

# plt.xlabel("Energy (eV)", fontsize = 12)
# plt.ylabel("Intensity (Arb. units)", fontsize = 12)
# plt.scatter(dataBetaPL[0],dataBetaPL[1]/max(dataBetaPL[1]), label = LabList[1] + ' 23K', color = 'darkcyan', s = 1)
# plt.plot(dataBetaPL[0],dataBetaPL[2]/max(dataBetaPL[1]), label = '3.5eV', color = GaussColours[6])
# plt.plot(dataBetaPL[0],dataBetaPL[3]/max(dataBetaPL[1]), label = '2.9eV', color = GaussColours[2])
# plt.plot(dataBetaPL[0],dataBetaPL[4]/max(dataBetaPL[1]), label = '2.5eV', color = GaussColours[1])
# plt.plot(dataBetaPL[0],dataBetaPL[5]/max(dataBetaPL[1]), color = 'k', linestyle = '--')
# plt.plot(dataBetaPL[0],dataBetaPL[6]/max(dataBetaPL[1]), color = 'k', linestyle = '--')

# plt.legend(loc = 'best', fontsize = 12)

# plt.figure(2,dpi = GraphResolution)
# plt.minorticks_on()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# plt.yscale('log')
# # plt.xlim(1.78622,4.5)
# plt.ylim(2e-3,1)

# plt.xlabel("Energy (eV)", fontsize = 12)
# plt.ylabel("Intensity (Arb. units)", fontsize = 12)
# plt.scatter(dataKappaPL[0],dataKappaPL[1]/max(dataKappaPL[1]), label = LabList[2] + ' 23K', color = 'darkcyan',s = 1)
# plt.plot(dataKappaPL[0],dataKappaPL[5]/max(dataKappaPL[1]), label = '2.1eV', color = 'orange')
# plt.plot(dataKappaPL[0],dataKappaPL[2]/max(dataKappaPL[1]), label = '2.6eV', color = GaussColours[2])
# plt.plot(dataKappaPL[0],dataKappaPL[3]/max(dataKappaPL[1]), label = '3.1eV', color = GaussColours[5])
# plt.plot(dataKappaPL[0],dataKappaPL[4]/max(dataKappaPL[1]), label = '3.6eV', color = GaussColours[6])
# plt.plot(dataKappaPL[0],dataKappaPL[6]/max(dataKappaPL[1]), color = 'k', linestyle = '--')
# plt.legend(loc = 'best', fontsize = 12)

# TickSize = 34
# LabelSize = 38
# LegSize = 30
# LineW = 3
# MarkSize = 8

# sizeFigSqaure = 10


TickSize = 28
LabelSize = 30
LegSize = 28
LineW = 2
MarkSize = 4

sizeFigSqaure = 30

plt.figure(6,dpi = GraphResolution, figsize = ((sizeFigSqaure + int(sizeFigSqaure/2)), sizeFigSqaure))
plt.minorticks_on()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.yscale('log')
plt.xlim(1.78622,4.5)

plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)

plt.ylim(2e-3,1)

plt.xlabel("Energy (eV)", fontsize =  LabelSize)
plt.ylabel("Intensity (Arb. units)", fontsize =  LabelSize)

plt.scatter(dataAlphaPLlowT[0],dataAlphaPLlowT[1]/max(dataAlphaPLlowT[1]), label = LabList[0] + ' 23K', color = 'darkcyan', s = MarkSize)
plt.plot(dataAlphaPLlowT[0],dataAlphaPLlowT[3]/max(dataAlphaPLlowT[1]), label = '3.9eV', color = GaussColours[6], lw = LineW)
plt.plot(dataAlphaPLlowT[0],dataAlphaPLlowT[4]/max(dataAlphaPLlowT[1]), label = '3.6eV', color = GaussColours[5], lw = LineW)
plt.plot(dataAlphaPLlowT[0],dataAlphaPLlowT[5]/max(dataAlphaPLlowT[1]), label = '3.1eV', color = GaussColours[3], lw = LineW)
plt.plot(dataAlphaPLlowT[0],dataAlphaPLlowT[6]/max(dataAlphaPLlowT[1]), label = '2.6eV', color = GaussColours[2], lw = LineW)
plt.plot(dataAlphaPLlowT[0],dataAlphaPLlowT[7]/max(dataAlphaPLlowT[1]), color = 'k', linestyle = '--', lw = LineW)
plt.legend(loc = 'best', fontsize = LegSize)

plt.figure(7,dpi = GraphResolution, figsize = (sizeFigSqaure * 1.5, sizeFigSqaure))
plt.minorticks_on()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.yscale('log')
plt.xlim(1.78622,4.5)
plt.ylim(2e-3,1)

plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)

plt.xlabel("Energy (eV)", fontsize =  LabelSize)
plt.ylabel("Intensity (Arb. units)", fontsize =  LabelSize)
plt.scatter(dataBetaPLlowT[0],dataBetaPLlowT[1]/max(dataBetaPLlowT[1]), label = LabList[1] + ' 23K', color = 'darkcyan', s = MarkSize)
plt.plot(dataBetaPLlowT[0],dataBetaPLlowT[3]/max(dataBetaPLlowT[1]), label = '3.6eV', color = GaussColours[6], lw = LineW)
plt.plot(dataBetaPLlowT[0],dataBetaPLlowT[4]/max(dataBetaPLlowT[1]), label = '3.3eV', color = GaussColours[5], lw = LineW)
plt.plot(dataBetaPLlowT[0],dataBetaPLlowT[5]/max(dataBetaPLlowT[1]), label = '2.9eV', color = GaussColours[3], lw = LineW)
plt.plot(dataBetaPLlowT[0],dataBetaPLlowT[6]/max(dataBetaPLlowT[1]), label = '2.5eV', color = GaussColours[2], lw = LineW)
plt.plot(dataBetaPLlowT[0],dataBetaPLlowT[7]/max(dataBetaPLlowT[1]), color = 'k', linestyle = '--', lw = LineW)

plt.legend(loc = 'best', fontsize = LegSize)

plt.figure(8,dpi = GraphResolution, figsize = (sizeFigSqaure * 1.5, sizeFigSqaure))
plt.minorticks_on()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.yscale('log')
plt.xlim(1.78622,4.5)
plt.ylim(2e-3,1)

plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)

plt.xlabel("Energy (eV)", fontsize =  LabelSize)
plt.ylabel("Intensity (Arb. units)", fontsize =  LabelSize)
plt.scatter(dataKappaPLlowT[0],dataKappaPLlowT[1]/max(dataKappaPLlowT[1]), label = LabList[2] + ' 23K', color = 'darkcyan',s = MarkSize)
plt.plot(dataKappaPLlowT[0],dataKappaPLlowT[3]/max(dataKappaPLlowT[1]), label = '2.2eV', color = 'orange', lw = LineW)
plt.plot(dataKappaPLlowT[0],dataKappaPLlowT[4]/max(dataKappaPLlowT[1]), label = '2.7eV', color = GaussColours[2], lw = LineW)
plt.plot(dataKappaPLlowT[0],dataKappaPLlowT[5]/max(dataKappaPLlowT[1]), label = '3.2eV', color = GaussColours[5], lw = LineW)
plt.plot(dataKappaPLlowT[0],dataKappaPLlowT[6]/max(dataKappaPLlowT[1]), label = '3.7eV', color = GaussColours[6], lw = LineW)
plt.plot(dataKappaPLlowT[0],dataKappaPLlowT[7]/max(dataKappaPLlowT[1]), color = 'k', linestyle = '--', lw = LineW)
plt.legend(loc = 'best', fontsize = LegSize)