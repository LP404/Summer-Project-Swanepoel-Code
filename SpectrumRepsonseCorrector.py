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

PhaseImport = 0

xMin = 280
xMax = 700

sizeFigSqaure = F1.CTI(20)

if PhaseImport == 0:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('SpectrumRepsonseCorrector.py')) + '\\PLTXT\\Hidden\\D'))
elif PhaseImport == 1:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('SpectrumRepsonseCorrector.py')) + '\\PLTXT\\Hidden\\1534'))
else:
    path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('SpectrumRepsonseCorrector.py')) + '\\PLTXT\\Hidden\\1198'))



CorrectionNames, CorrectionData = nPF.txtGrabber('SpectrumRepsonseCorrector.py', '\\SystemResponse\\SystemResponseScripts\\PL\\ResponseTxt\\filt', False)

CorrectionNamesUnfilt, CorrectionDataUnfilt = nPF.txtGrabber('SpectrumRepsonseCorrector.py', '\\SystemResponse\\SystemResponseScripts\\PL\\ResponseTxt\\unfilt', False)

AvgCor = CorrectionData[0]
LongCor = CorrectionData[1]

AvgCorTrim = F.Trim(AvgCor,xMin,xMax)
LongCorTrim = F.Trim(LongCor,xMin,xMax)

# AvgCorTrim_eV = ConvertDependance(AvgCorTrim[0],AvgCorTrim[1])
# LongCorTrim_eV = ConvertDependance(LongCorTrim[0],LongCorTrim[1])

# AvgCorTrim_eV = np.flip(AvgCorTrim_eV,1)

# LongCorTrim_eV = np.flip(LongCorTrim_eV,1)



AvgCorUnfilt = CorrectionDataUnfilt[0]
LongCorUnfilt = CorrectionDataUnfilt[1]

AvgCorTrimUnfilt = F.Trim(AvgCorUnfilt,xMin,xMax)
LongCorTrimUnfilt = F.Trim(LongCorUnfilt,xMin,xMax)

# AvgCorTrim_eVUnfilt = ConvertDependance(AvgCorTrimUnfilt[0],AvgCorTrimUnfilt[1])
# LongCorTrim_eVUnfilt = ConvertDependance(LongCorTrimUnfilt[0],LongCorTrimUnfilt[1])

# AvgCorTrim_eVUnfilt = np.flip(AvgCorTrim_eVUnfilt,1)

# LongCorTrim_eVUnfilt = np.flip(LongCorTrim_eVUnfilt,1)


files = natsorted(files)


#In nm

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i][:-4]
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",", skiprows=2).T
    #, skiprows=2
    vars()[files[i]][1] = vars()[files[i]][1]

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)


    #Change 0.075 back to 0.2
    # vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])
    vars()['yFiltered'+str(i)] = vars()[files[i]+'T'][1]


Red = np.linspace(0,255,len(files)+1).astype('int32')
Blue = np.flip(Red)


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple']
LineStyle = ['solid','dashed']

TickSize = 14
LabelSize = 14
LegSize = 10




# plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
plt.figure(0,figsize=(1.5*CTI(16),CTI(16)), dpi=600)
for i in range(len(files)):
# for i in range(0,1):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.minorticks_on()
#    plt.grid(which='major', color='k', linestyle='-')
#    plt.grid(which='minor', color='darkgray', linestyle='--')
    
    xVals = vars()[files[i]+'T'][0]
    yValsUncor = vars()['yFiltered'+str(i)]

    # xVals,yVals = ConvertDependance(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])

    yValCorAvg = np.interp(xVals,AvgCorTrim[0],AvgCorTrim[1])
    yValCorLong = np.interp(xVals,LongCorTrim[0],LongCorTrim[1])
    
    yValCorAvgUnfilt = np.interp(xVals,AvgCorTrimUnfilt[0],AvgCorTrimUnfilt[1])
    yValCorLongUnfilt = np.interp(xVals,LongCorTrimUnfilt[0],LongCorTrimUnfilt[1])

    
    # xVals = np.flip(xVals)
    # yVals = np.flip(yVals)
    

    # yValCorAvg = np.interp(xVals,AvgCorTrim_eV[0],AvgCorTrim_eV[1])
    # yValCorLong = np.interp(xVals,LongCorTrim_eV[0],LongCorTrim_eV[1])
    
    # yVals = yValsUncor
    # yVals2 = yValsUncor
    
    yVals = yValsUncor/yValCorAvg
    yVals2 = yValsUncor/yValCorAvgUnfilt

    # yVals = yValsUncor/yValCorLong
    # yVals = yValsUncor/yValCorLongUnfilt
    
    xVals_eV, yVals_eV = ConvertDependance(xVals,yVals)
    xVals2_eV, yVals2_eV = ConvertDependance(xVals,yVals2)
    
    xVals_eV = np.flip(xVals_eV)
    yVals_eV = np.flip(yVals_eV)
    
    xVals2_eV = np.flip(xVals2_eV)
    yVals2_eV = np.flip(yVals2_eV)
    
    

    
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

    # plt.semilogy(xVals,yVals, label = files[i], color = (Red[i]/255,0,Blue[i]/255))
    
    # plt.semilogy(xVals,yVals2, label = files[i], color = (Red[i]/255,0,Blue[i]/255))

    
    plt.semilogy(xVals_eV,yVals_eV, label = files[i], color = (Red[i]/255,0,Blue[i]/255))
    
    # plt.semilogy(xVals_eV,yVals2_eV, label = files[i], color = (Red[i]/255,0,Blue[i]/255))

    # if PhaseImport == 0:
    #     np.savetxt(f'PLTXT\\txtOut\\b\\β-Ga₂O₃ {files[i]} Specrtra.txt',np.c_[xVals_eV,yVals_eV], delimiter=",", fmt='%.3f')
    #     np.savetxt(f'PLTXT\\txtOut\\b\\β-Ga₂O₃ {files[i]} Specrtra Unfiltered.txt',np.c_[xVals2_eV,yVals2_eV], delimiter=",", fmt='%.3f')
    # elif PhaseImport == 1:
    #     np.savetxt(f'PLTXT\\txtOut\\a\\α-Ga₂O₃ {files[i]} Specrtra.txt',np.c_[xVals_eV,yVals_eV], delimiter=",", fmt='%.3f')
    #     np.savetxt(f'PLTXT\\txtOut\\a\\α-Ga₂O₃ {files[i]} Specrtra Unfiltered.txt',np.c_[xVals2_eV,yVals2_eV], delimiter=",", fmt='%.3f')
    # else:
    #     np.savetxt(f'PLTXT\\txtOut\\k\\κ-Ga₂O₃ {files[i]} Specrtra.txt',np.c_[xVals_eV,yVals_eV], delimiter=",", fmt='%.3f')
    #     np.savetxt(f'PLTXT\\txtOut\\k\\κ-Ga₂O₃ {files[i]} Specrtra Unfiltered.txt',np.c_[xVals2_eV,yVals2_eV], delimiter=",", fmt='%.3f')




if PhaseImport == 0:
    plt.title('Photoluminescence of β-Ga₂O₃', fontsize = LabelSize)
    # plt.axvline(x = 2.5, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
    # plt.axvline(x = 3.2, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)
    # plt.axvline(x = 3.4, ymin = 0, ymax = 1 , color = "purple", linewidth = 3)
    # plt.axvline(x = 3.8, ymin = 0, ymax = 1 , color = "darkviolet", linewidth = 3)
elif PhaseImport == 1:
    plt.title('Photoluminescence of α-Ga₂O₃', fontsize = LabelSize) 
    # plt.axvline(x = 2.5, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
    # plt.axvline(x = 2.9, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)
    # plt.axvline(x = 3.3, ymin = 0, ymax = 1 , color = "purple", linewidth = 3)
    # plt.axvline(x = 3.7, ymin = 0, ymax = 1 , color = "darkviolet", linewidth = 3)
else:
    plt.title('Photoluminescence of κ-Ga₂O₃', fontsize = LabelSize) 
    # plt.axvline(x = 2.0, ymin = 0, ymax = 1 , color = "crimson", linewidth = 3)
    # plt.axvline(x = 2.4, ymin = 0, ymax = 1 , color = "green", linewidth = 3)
    # plt.axvline(x = 2.9, ymin = 0, ymax = 1 , color = "deepskyblue", linewidth = 3)    
    

plt.xlabel("Energy (eV)", fontsize = LabelSize)
# plt.xlabel("Wavelength (nm)", fontsize = LabelSize)
plt.ylabel("Intensity (A.U.)", fontsize = LabelSize)
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
plt.xlim(1.77,4.43)
# plt.xlim(280,700)
# plt.ylim(yMinVals,yMaxVals)
plt.legend(fontsize = LegSize, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)