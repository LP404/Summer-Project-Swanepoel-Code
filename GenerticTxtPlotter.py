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

PhaseImport = 2

Slitwidth = 2

xMin = 280
xMax = 700

sizeFigSqaure = F1.CTI(20)

# path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\Data'))

# path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\FitOptTest\\raw'))


if PhaseImport == 0:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Hidden\\D'))
elif PhaseImport == 1:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Hidden\\1534'))
else:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Hidden\\1198'))

# path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Poster\\α'))
# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Poster\\β'))
# path2, dirs2, files2 = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\PLTXT\\Poster\\κ'))

# headerAlphaCL, dataAlphaCL = nPF.datGrabber('ExportedDataAlpha','fitykExports\\')
# headerBetaCL, dataBetaCL = nPF.datGrabber('ExportedDataBeta','fitykExports\\')
# headerKappaCL, dataKappaCL = nPF.datGrabber('ExportedDataKappa','fitykExports\\')

# headerAlphaCLLT, dataAlphaCLLT = nPF.datGrabber('MeanSpectra10keV','fitykExports\\')
# TheaderBetaCLLT, dataBetaCLLT = nPF.datGrabber('BetaMeanSpectra','fitykExports\\')
# headerKappaCLLT, dataKappaCLLT = nPF.datGrabber('MeanKappaSpecrtra','fitykExports\\')

CorrectionNames, CorrectionData = nPF.txtGrabber('GenerticTxtPlotter.py', '\\SystemResponse\\SystemResponseScripts\\PL\\ResponseTxt\\filt', False)

CorrectionNamesUnfilt, CorrectionDataUnfilt = nPF.txtGrabber('GenerticTxtPlotter.py', '\\SystemResponse\\SystemResponseScripts\\PL\\ResponseTxt\\unfilt', False)

AvgCor = CorrectionData[0]
LongCor = CorrectionData[1]

AvgCorTrim = F.Trim(AvgCor,xMin,xMax)
LongCorTrim = F.Trim(LongCor,xMin,xMax)


# AvgCorTrim_eV = ConvertDependance(AvgCorTrim[0],AvgCorTrim[1])
# LongCorTrim_eV = ConvertDependance(LongCorTrim[0],LongCorTrim[1])

# AvgCorTrim_eV = np.flip(AvgCorTrim_eV,1)

# LongCorTrim_eV = np.flip(LongCorTrim_eV,1)


files = natsorted(files)

# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotter.py')) + '\\Data\\RetakeWithPaint\\Mugove\\Substrate'))

# path2, dirs2, files2, data2, header2 = nPF.fileGrabber('GenerticTxtPlotter.py','\\FitOptTest\\csv','.csv')

#In nm

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i][:-4]
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",", skiprows=2).T
    #, skiprows=2
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)


    #Change 0.075 back to 0.2
    vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])



    # files1[i] = files1[i][:-4]
    
    # vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",").T
    # vars()[files1[i]][1] = vars()[files1[i]][1] / 100

    # #T stands for Truncated
    # vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],xMin,xMax)


    # #Change 0.075 back to 0.2
    # vars()['yFiltered1'+str(i)] = F.NoiseFilter(3,0.2,vars()[files1[i]+'T'][1])

    # files2[i] = files2[i][:-4]
    
    # vars()[files2[i]] = np.loadtxt(open(path2 + "\\" + files2[i] + ".txt", "rb"), delimiter=",").T
    # vars()[files2[i]][1] = vars()[files2[i]][1] / 100

    # #T stands for Truncated
    # vars()[files2[i]+'T'] = F.Trim(vars()[files2[i]],xMin,xMax)


    # #Change 0.075 back to 0.2
    # vars()['yFiltered2'+str(i)] = F.NoiseFilter(3,0.2,vars()[files2[i]+'T'][1])    
    
    # if i == 0:
    #     #All the xValues are the same, only need to do this once
    #     xP = np.linspace(min(vars()[files[i]+'T'][0]),max(vars()[files[i]+'T'][0]),10001)
    # else:
    #     pass



# for i in range(len(files1)):
#     files1[i] = files1[i][:-4]
    
#     vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    
#     #This can be commented out with an alternet approach
#     vars()[files1[i]][0] = vars()[files1[i]][0]
    
#     vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
#     #vars()[files1[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files1[i]],xMin,xMax))
#     vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],xMin,xMax)

#     if i == 0:
#         #All the xValues are the same, only need to do this once
#         xP = np.linspace(min(vars()[files1[i]+'T'][0]),max(vars()[files1[i]+'T'][0]),10001)
#     else:
#         pass

#     vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])


Red = np.linspace(0,255,len(files)+1).astype('int32')
Blue = np.flip(Red)


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple']
LineStyle = ['solid','dashed']


# plt.figure(0,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# plt.title('Transmittance of κ-Ga₂O₃',fontsize = 18)
# plt.xlabel("Wavelength (nm)",fontsize = 18)
# plt.ylabel("Transmission (%)",fontsize = 18)
# plt.minorticks_on()
# plt.xticks(fontsize = 18)
# plt.yticks(fontsize = 18)
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')    
#for i in range(len(files)):
# plt.plot(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1]*100, label = LabList[i], color = LineColour[i], linestyle = LineStyle[0])

#plt.legend(fontsize = 18)


TickSize = 14
LabelSize = 14
LegSize = 10




# plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
plt.figure(0,figsize=(1.5*CTI(16),CTI(16)), dpi=600)
for i in range(len(files)):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.minorticks_on()
#    plt.grid(which='major', color='k', linestyle='-')
#    plt.grid(which='minor', color='darkgray', linestyle='--')
    
    xVals = vars()[files[i]+'T'][0]
    yVals = vars()['yFiltered'+str(i)] * 100

    yValCorAvg = np.interp(xVals,AvgCorTrim[0],AvgCorTrim[1])
    yValCorLong = np.interp(xVals,LongCorTrim[0],LongCorTrim[1])


    # xVals,yVals = ConvertDependance(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)]*100)
    
    
    # xVals = np.flip(xVals)
    # yVals = np.flip(yVals)
    

    # yValCorAvg = np.interp(xVals,AvgCorTrim_eV[0],AvgCorTrim_eV[1])
    # yValCorLong = np.interp(xVals,LongCorTrim_eV[0],LongCorTrim_eV[1])
    

    
    yVals = yVals/yValCorAvg
    # yVals = yVals/yValCorLong
    
    if PhaseImport == 0:
        if i == 0:
            yMaxVals = max(yVals) * 1.25
        elif i == (len(files) - 1):
            yMinVals = min(yVals) * 5
        else:
            pass  
    elif PhaseImport == 1:
        if i == 0:
            yMaxVals = max(yVals) * 1.25
        elif i == (len(files) - 1):
            yMinVals = min(yVals) * 5
        else:
            pass
    else:
        if i == 0:
            yMaxVals = max(yVals) * 1.25
        elif i == (len(files) - 1):
            yMinVals = min(yVals) * 4
        else:
            pass

    
    plt.semilogy(xVals,yVals, label = files[i], color = (Red[i]/255,0,Blue[i]/255))
    # np.savetxt('CorrectedSpecrtra.txt',np.c_[SaveNewX,SaveNewY], delimiter=",", fmt='%.3f')
    # np.savetxt('ResponseOutputLong.txt',np.c_[SaveNewX,SaveNewY2], delimiter=",", fmt='%.3f')
    # plt.semilogy(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)], label = LabList[0] + ' ' + files[i], color = LineColour[0],linestyle = LineStyle[i])
    # plt.semilogy(vars()[files1[i]+'T'][0],vars()['yFiltered1'+str(i)], label = LabList[1] + ' ' + files1[i], color = LineColour[1],linestyle = LineStyle[i])
    # plt.semilogy(vars()[files2[i]+'T'][0],vars()['yFiltered2'+str(i)], label = LabList[2] + ' ' + files2[i], color = LineColour[2],linestyle = LineStyle[i])
    # # plt.title(f'{files[i]}')
## ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']

if PhaseImport == 0:
    plt.title('Photoluminescence of β-Ga₂O₃', fontsize = LabelSize)
    plt.axvline(x = 2.5, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
    plt.axvline(x = 3.2, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)
    plt.axvline(x = 3.4, ymin = 0, ymax = 1 , color = "purple", linewidth = 3)
    plt.axvline(x = 3.8, ymin = 0, ymax = 1 , color = "darkviolet", linewidth = 3)
elif PhaseImport == 1:
   plt.title('Photoluminescence of α-Ga₂O₃', fontsize = LabelSize) 
   plt.axvline(x = 2.5, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
   plt.axvline(x = 2.9, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)
   plt.axvline(x = 3.3, ymin = 0, ymax = 1 , color = "purple", linewidth = 3)
   plt.axvline(x = 3.7, ymin = 0, ymax = 1 , color = "darkviolet", linewidth = 3)
else:
   plt.title('Photoluminescence of κ-Ga₂O₃', fontsize = LabelSize) 
   plt.axvline(x = 2.0, ymin = 0, ymax = 1 , color = "crimson", linewidth = 3)
   plt.axvline(x = 2.4, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
   plt.axvline(x = 2.9, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)    
    

plt.xlabel("Energy (eV)", fontsize = LabelSize)
plt.ylabel("Intensity (A.U.)", fontsize = LabelSize)
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
# plt.xlim(1.77,4.43)
plt.xlim(280,700)
# plt.ylim(yMinVals,yMaxVals)
plt.legend(fontsize = LegSize, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)

# plt.figure(1,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# for i in range(len(files)):
#     vars()['yFiltered Energy'+str(i)] = vars()['yFiltered'+str(i)] * ((vars()[files[i]+'T'][0]**2) / (6.62607015e-34 * 299792458))
#     vars()['energyXscale'+str(i)] = ((6.62607015e-34 * 299792458) / (vars()[files[i]+'T'][0]*1e-9)) * 6.242e18
    
#     vars()['yFiltered Energy1'+str(i)] = vars()['yFiltered1'+str(i)] * ((vars()[files1[i]+'T'][0]**2) / (6.62607015e-34 * 299792458))
#     vars()['energyXscale1'+str(i)] = ((6.62607015e-34 * 299792458) / (vars()[files1[i]+'T'][0]*1e-9)) * 6.242e18
    
#     vars()['yFiltered Energy2'+str(i)] = vars()['yFiltered2'+str(i)] * ((vars()[files2[i]+'T'][0]**2) / (6.62607015e-34 * 299792458))
#     vars()['energyXscale2'+str(i)] = ((6.62607015e-34 * 299792458) / (vars()[files2[i]+'T'][0]*1e-9)) * 6.242e18



#     plt.minorticks_on()
#     plt.grid(which='major', color='k', linestyle='-')
#     plt.grid(which='minor', color='darkgray', linestyle='--')
#     plt.semilogy(vars()['energyXscale'+str(i)],vars()['yFiltered Energy'+str(i)], label = LabList[0] + ' ' + files[i], color = LineColour[0],linestyle = LineStyle[i])
#     plt.semilogy(vars()['energyXscale1'+str(i)],vars()['yFiltered Energy1'+str(i)], label = LabList[1] + ' ' + files1[i], color = LineColour[1],linestyle = LineStyle[i])
#     plt.semilogy(vars()['energyXscale2'+str(i)],vars()['yFiltered Energy2'+str(i)], label = LabList[2] + ' ' + files2[i], color = LineColour[2],linestyle = LineStyle[i])

    
#     plt.title(f'Photoluminescence of Ga₂O₃ polymorphs', fontsize = 12)
#     plt.xlabel("Energy (eV)", fontsize = 12)
#     plt.ylabel("Intensity (counts/eV)", fontsize = 12)
#     plt.legend(loc = 'best', fontsize = 8)
    
# GraphResolution = 100
# sizeFigSqaure = CTI(30)

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

# plt.figure(1,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# for i in range(len(files1)):
#     # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
#     plt.minorticks_on()
#     plt.grid(which='major', color='k', linestyle='-')
#     plt.grid(which='minor', color='darkgray', linestyle='--')
#     plt.plot(xP,vars()['yP'+files1[i]], label = files1[i])
#     # plt.title(f'{files[i]}')
#     plt.title(f'Substrate Transmittnce Uncorrected')
#     plt.xlabel("Wavelength (nm)")
#     plt.ylabel("Transmission/Absorbance")
#     plt.legend()
