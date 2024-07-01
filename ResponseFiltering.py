import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import nPlotterFunctions as nPF
import os
import SpectralResponseInterpArray as SRIA
import time
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi, savgol_filter
from math import log10, floor
import copy
import random
from scipy.ndimage import gaussian_filter1d
from scipy import interpolate
from whittaker_eilers import WhittakerSmoother
from sklearn.metrics import r2_score

def TrimSingle(Array,start,end):
    delPoints = np.where((Array < start) | (Array > end))[0]
    
    newArray = np.delete(Array,delPoints)
    
    return newArray


def Binning(array, n):
    return array[:(array.size // n) * n].reshape(-1, n).mean(axis=1)
    
def Man(Jonkles):
    ExpArr = np.floor(np.log10(abs(Jonkles)))
    
    Loc = np.where(ExpArr > (np.mean(ExpArr) + 3))[0]
    
    if Loc.size == 0:
        pass
    
    else:
        
        if Loc[0] == 0:
            Val = ExpArr[0] - ExpArr[1]
            Div = 10**Val
            Jonkles[0] = Jonkles[0] / Div
        
        if Loc[-1] == (len(ExpArr)-1):
            Val = ExpArr[-1] - ExpArr[-2]
            Div = 10**Val
            Jonkles[-1] = Jonkles[-1] / Div
    
        for i in range(len(Loc)):
            Val = ExpArr[Loc[i]] - np.floor(np.mean(ExpArr))
            Div = 10**Val
            Jonkles[Loc[i]] = Jonkles[Loc[i]] / Div
    
    
    return Jonkles

def RemoveZero(Array):
    
    if Array[0] == 0:
        Array[0] == np.where(Array != 0)[0][0]
    if Array[-1] == 0:
        Array[-1] == np.where(Array != 0)[0][-1]
    
    Locs = np.where(Array == 0)[0]
    for i in range(len(Locs)):
        Array[Locs[i]] = ((Array[(Locs[i]+1)] + Array[(Locs[i]-1)]) /2)
        
    return Array
        


def NoiseFilter(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    yFiltered = filtfilt(b, a, yArray)
    
    return yFiltered

def NoiseFilter2(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    zi = lfilter_zi(b, a)
    z, _ = lfilter(b, a, yArray, zi=zi*yArray[0])
    z2, _ = lfilter(b, a, z, zi=zi*z[0])
    
    return z2

def NoiseFilter3(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    zi = lfilter_zi(b, a)
    z, _ = lfilter(b, a, yArray, zi=zi*yArray[0])
    z2, _ = lfilter(b, a, z, zi=zi*z[0])
    yFiltered = filtfilt(b, a, z2)
    
    return yFiltered

def noise_filter(order, cutoff_frequency, signal, initial_conditions=None, filter_type='lowpass'):
    """
    Apply Butterworth filter to remove noise from a signal.

    Parameters:
    - order (int): Filter order.
    - cutoff_frequency (float): Cutoff frequency of the filter.
    - signal (array-like): Input signal.
    - initial_conditions (array-like, optional): Initial conditions for the filter.
    - filter_type (str, optional): Type of filter ('lowpass', 'highpass', 'bandpass', etc.).

    Returns:
    - array-like: Filtered signal.
    """

    # Validate parameters
    if not isinstance(order, int) or order <= 0:
        raise ValueError("Order must be a positive integer.")
    if not isinstance(cutoff_frequency, (int, float)) or cutoff_frequency <= 0:
        raise ValueError("Cutoff frequency must be a positive number.")

    # Design Butterworth filter
    b, a = butter(order, cutoff_frequency, btype=filter_type)

    # Calculate initial conditions if not provided
    if initial_conditions is None:
        initial_conditions = lfilter_zi(b, a)

    # Apply first lfilter
    filtered_signal_1, _ = lfilter(b, a, signal, zi=initial_conditions)

    # Apply second lfilter based on the result of the first one
    filtered_signal_2, _ = lfilter(b, a, filtered_signal_1, zi=initial_conditions)

    return filtered_signal_2


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


def SpikeCorrect(xArray,yArray,xStart,xEnd,samples):
    Locs1 = np.where((xArray >= (xStart - samples)) & (xArray <= xStart))[0]
    Locs2 = np.where((xArray >= xEnd) & (xArray <= (xEnd + samples)))[0]
    Locs1 = np.append(Locs1,Locs2)  

    meanY = np.mean(yArray[Locs1])
    maxY = max(yArray[Locs1])
    minY = min(yArray[Locs1])
    
    Locs3 = np.where((xArray >= xStart) & (xArray <= xEnd))[0]
    TargetY = yArray[Locs3]
    
    for i in range(len(TargetY)):
        if TargetY[i] <= maxY and TargetY[i] >= minY:
            pass
        else:
            EndCon = 0
            while EndCon == 0:
                if TargetY[i] > maxY:
                    div = TargetY[i] / meanY
                    divNew = div + random.uniform(-0.25,0.25)
                    TargetY[i] = meanY * divNew
                elif TargetY[i] <  minY:
                    div = TargetY[i] / meanY
                    divNew = div + random.uniform(-0.25,0.25)
                    TargetY[i] = meanY * divNew                    
                else:
                    EndCon = 1
    
    newYarr = copy.deepcopy(yArray)
    
    for i in range(len(TargetY)):
        newYarr[Locs3[i]] = TargetY[i]
                
    
    return newYarr
    
def Rsquared(y,yFit):
    
    rVal = 1
    
    return rVal

CLFileNames, CLData = nPF.txtGrabber('SpectralResponse.py', '\\SystemResponse\\SystemResponseScripts\\CL', False)


pathQuanta1,dirsQuanta1,filesQuanta1,dataQuanta1  = nPF.fileGrabberCSVNH('SpectralResponse.py', "\\SystemResponse\\SystemResponseScripts\\PL\\LampInQuanta\\OceanOpticsLamp")
pathQuantaAvg,dirsQuantaAvg,filesQuantaAvg,dataQuantaAvg  = nPF.fileGrabberCSVNH('SpectralResponse.py', "\\SystemResponse\\SystemResponseScripts\\PL\\LampInQuanta\\OceanOpticsLamp\\Avg")
AndorFileNames, AndorData = nPF.txtGrabber('SpectralResponse.py', '\\SystemResponse\\Andor', False)


diff = 0
diffArray = np.array([])
for i in range(1,len(CLData[0][0])):
    diff += abs((CLData[0][0][i] - CLData[0][0][i-1]))
    diffArray = np.append(diffArray,abs((CLData[0][0][i] - CLData[0][0][i-1])))
    

AvgDiff = diff / len(CLData[0][0])

dt = AvgDiff

y = CLData[0][1]

yFilt = NoiseFilter(3,0.1,y)
yFilt2 = NoiseFilter2(3,0.1,y)
yFilt3 = NoiseFilter3(3,0.1,y)

yFilt4 = gaussian_filter1d(y,1)
yFilt5 = gaussian_filter1d(y,3)
yFilt6 = gaussian_filter1d(y,6)
yFilt7 = gaussian_filter1d(y,10)
yFilt8 = gaussian_filter1d(y,12)

n = len(y)
fhat = np.fft.fft(y,n)

PSD = np.zeros_like(fhat)

for i in range(len(PSD)):
    A = fhat[i] * np.conj(fhat[i])
    B = A/n
    PSD[i] = B

freq = (1/(dt*n)) * np.arange(n)
L = np.arange(1,np.floor(n/2),dtype = 'int')

fig,ax = plt.subplots(2,1)

fig.tight_layout()

plt.sca(ax[0])
plt.plot(CLData[0][0],CLData[0][1], linewidth=1.5, label = 'Noisy')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])

plt.sca(ax[1])
plt.plot(freq[L],PSD[L], linewidth=2, label = 'Noisy')
plt.xlim(0,0.02)


inds = PSD > 1e20
PSDcleaned = PSD * inds
fhat = inds * fhat
ffilt = np.fft.ifft(fhat)




fig1,ax1 = plt.subplots(11,1,figsize=(25,35))
fig1.tight_layout()

plt.sca(ax1[0])
plt.plot(CLData[0][0],CLData[0][1], linewidth=1, label = 'Noisy')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[1])
plt.plot(CLData[0][0],ffilt, linewidth=1, label = 'Cleaned by FFT')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[2])
plt.plot(CLData[0][0],yFilt, linewidth=1, label = 'Cleaned by filtfilt')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[3])
plt.plot(CLData[0][0],yFilt2, linewidth=1, label = 'Cleaned by lfilter twice')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()


plt.sca(ax1[4])
plt.plot(CLData[0][0],yFilt3, linewidth=1, label = 'Cleaned by both')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()


plt.sca(ax1[5])
plt.plot(CLData[0][0],yFilt4, linewidth=1, label = 'Gaussian Filter sigma = 1')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[6])
plt.plot(CLData[0][0],yFilt5, linewidth=1, label = 'Gaussian Filter sigma = 3')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[7])
plt.plot(CLData[0][0],yFilt6, linewidth=1, label = 'Gaussian Filter sigma = 6')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[8])
plt.plot(CLData[0][0],yFilt7, linewidth=1, label = 'Gaussian Filter sigma = 10')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[9])
plt.plot(CLData[0][0],yFilt8, linewidth=1, label = 'Gaussian Filter sigma = 15')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()

plt.sca(ax1[10])
plt.scatter(CLData[0][0],CLData[0][1], linewidth=1, label = 'Noisy')
plt.plot(CLData[0][0],yFilt8, linewidth=1, label = 'Gaussian Filter sigma = 15')
plt.xlim(CLData[0][0][0],CLData[0][0][-1])
plt.legend()


plt.show()


for i in range(len(dataQuantaAvg)):
    vars()[f"DataSet {[i]} x"] = np.flip(np.array(nPF.ListExtract(dataQuantaAvg[i],0)))
    vars()[f"DataSet {[i]} y"] = np.flip(np.array(nPF.ListExtract(dataQuantaAvg[i],1)))
    

AvgDataSetY = np.zeros_like(vars()[f"DataSet {[0]} y"])

for i in range(len(vars()[f"DataSet {[0]} x"])):
    total = 0
    for j in range(len(dataQuantaAvg)):
        total += vars()[f"DataSet {[j]} y"][i]
        mean = (total / len(dataQuantaAvg))
        
    AvgDataSetY[i] = mean


LongTakeQuantaX = np.flip(np.array(nPF.ListExtract(dataQuanta1[0],0)))
LongTakeQuantaY = np.flip(np.array(nPF.ListExtract(dataQuanta1[0],1)))

AvgDataSetY = RemoveZero(AvgDataSetY)
LongTakeQuantaY = RemoveZero(LongTakeQuantaY)


NewYforXSysResp = np.interp(vars()[f"DataSet {[0]} x"],CLData[0][0],CLData[0][1])

LongTakeQuantaY = LongTakeQuantaY / max(LongTakeQuantaY)
AvgDataSetY = AvgDataSetY / max(AvgDataSetY)
NewYforXSysResp = NewYforXSysResp / max(NewYforXSysResp)

CorrectedData = AvgDataSetY / NewYforXSysResp

CorrectedData2 = LongTakeQuantaY / NewYforXSysResp

PLsetupX = AndorData[0][0]
PLsetupY = AndorData[0][1]

PLsetUPyNew = np.interp(vars()[f"DataSet {[0]} x"],AndorData[0][0],AndorData[0][1])

plt.plot(AndorData[0][0],AndorData[0][1])
plt.plot(vars()[f"DataSet {[0]} x"],PLsetUPyNew)

SystemResponse = PLsetUPyNew / CorrectedData
SystemResponse2 = PLsetUPyNew / CorrectedData2

SystemResponse = Man(SystemResponse)
SystemResponse2 = Man(SystemResponse2)

SystemResponseOriginal = copy.deepcopy(SystemResponse)
SystemResponseOriginal2 = copy.deepcopy(SystemResponse2)

plt.figure(6)
plt.plot(vars()[f"DataSet {[0]} x"],CorrectedData)
plt.plot(vars()[f"DataSet {[0]} x"],CorrectedData2)
plt.ylim(0,1)
plt.xlim(280,max(vars()[f"DataSet {[0]} x"]))

TickSize = 16
LabelSize = 16
LegSize = 14
BinCount = 3
ValA = 3
ValB = 0.15


SystemResponse = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse,650,665,10)
SystemResponse2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse2,650,665,10)

SystemResponse = Binning(SystemResponse,BinCount)
SystemResponse2 = Binning(SystemResponse2,BinCount)

SystemResponse = NoiseFilter(ValA,ValB,SystemResponse)
SystemResponse2 = NoiseFilter(ValA,ValB,SystemResponse2)

xBinned = Binning(vars()[f"DataSet {[0]} x"],BinCount)


# SystemResponseOriginal = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponseOriginal,650,665,10)
# SystemResponseOriginal2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponseOriginal2,650,665,10)

# SaveNewX = vars()[f"DataSet {[0]} x"][np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]
# SaveNewY = np.interp(SaveNewX,xBinned,SystemResponse)
# SaveNewY2 = np.interp(SaveNewX,xBinned,SystemResponse2)

# SaveNewYUnfil = np.interp(SaveNewX,vars()[f"DataSet {[0]} x"],SystemResponseOriginal)
# SaveNewY2Unfil = np.interp(SaveNewX,vars()[f"DataSet {[0]} x"],SystemResponseOriginal2)

# np.savetxt('ResponseOutputAvg.txt',np.c_[SaveNewX,SaveNewY], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLong.txt',np.c_[SaveNewX,SaveNewY2], delimiter=",", fmt='%.3f')

# np.savetxt('ResponseOutputAvgUnfilt.txt',np.c_[SaveNewX,SaveNewYUnfil], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLongUnfilt.txt',np.c_[SaveNewX,SaveNewY2Unfil], delimiter=",", fmt='%.3f')