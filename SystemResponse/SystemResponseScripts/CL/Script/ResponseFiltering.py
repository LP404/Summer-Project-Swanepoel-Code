import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi, savgol_filter
from scipy.ndimage import gaussian_filter1d
from whittaker_eilers import WhittakerSmoother
from sklearn.metrics import r2_score

def TrimSingle(Array,start,end):
    delPoints = np.where((Array < start) | (Array > end))[0]
    
    newArray = np.delete(Array,delPoints)
    
    return newArray


def Binning(array, n):
    return array[:(array.size // n) * n].reshape(-1, n).mean(axis=1)
    

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

             
def txtGrabber(scriptName,filePath,skipheader):
    
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))
    
    data = []
    
    #The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
    for i in range(len(fileNames)):
        fileNames[i] = fileNames[i][:-4]
        
        if skipheader == True:
            vars()[fileNames[i]] = np.loadtxt(open(path + "\\" + fileNames[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
        else:
            vars()[fileNames[i]] = np.loadtxt(open(path + "\\" + fileNames[i] + ".txt", "rb"), delimiter=",").T
            
        data.append(vars()[fileNames[i]])
    
    return fileNames, data

CLFileNames, CLData = txtGrabber('ResponseFiltering.py', '\\Input', False)

TickSize = 16
LabelSize = 16
LegSize = 14
BinCount = 3
ValA = 3
ValB = 0.15

STD = 12

window = 179
poly = 5

# lam = 100
# lam = 316227.7660168379
lam = 5e6
# lam2 = 3162.2776601683795

smoothOrder = 3
# smoothOrder = 2

yFilt9 = gaussian_filter1d(CLData[0][1],STD)
yFilt10 = savgol_filter(CLData[0][1],window,poly)

ws = WhittakerSmoother(lmbda=lam, order=smoothOrder, data_length = len(CLData[0][1]))

yFilt11 = ws.smooth(CLData[0][1])

yFilt12 = ws.smooth_optimal(CLData[0][1], break_serial_correlation=False)

optimally_smoothed_series = yFilt12.get_optimal().get_smoothed()
optimal_lambda = yFilt12.get_optimal().get_lambda()

plt.figure(8, figsize = (24,10))
plt.minorticks_on()
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
plt.scatter(CLData[0][0],CLData[0][1], linewidth=1, label = f'Raw Data, R^2 = {np.around(r2_score(CLData[0][1],CLData[0][1]),5)}', color = 'darkgrey')
plt.plot(CLData[0][0],yFilt9, linewidth=3, label = f'Gaussian Filter sigma = {STD},  R^2 = {np.around(r2_score(CLData[0][1],yFilt9),5)}', color = 'k', linestyle = '--')
plt.plot(CLData[0][0],yFilt10, linewidth=3, label = f'Savitzky-Golay Filter window = {window}, poly order = {poly},  R^2 = {np.around(r2_score(CLData[0][1],yFilt10),5)}', color = '#1a80bb', linestyle = '--')
plt.plot(CLData[0][0],yFilt11, linewidth=3, label = f'Whittaker-Eilers Filter Roughness penalty = {lam},\nSmoothing Order = {smoothOrder},  R^2 = {np.around(r2_score(CLData[0][1],yFilt11),5)}', color = '#ea801c', linestyle = '--')
plt.title('System response suggestion CL', fontsize = LabelSize)
plt.ylabel('System Response /nm (Arb units)', fontsize = LabelSize)
plt.xlabel("Wavelength (nm)", fontsize = LabelSize)
plt.legend(fontsize = LegSize)
plt.xlim(280,700)

plt.figure(9, figsize = (24,10))
plt.minorticks_on()
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
plt.scatter(CLData[0][0],CLData[0][1], linewidth=1, label = 'Raw Data', color = 'darkgrey')
plt.plot(CLData[0][0],yFilt9, linewidth=3, label = f'Gaussian Filter sigma = {STD},  R^2 = {np.around(r2_score(CLData[0][1],yFilt9),5)}', color = 'k', linestyle = '--')
plt.title('System response suggestion CL', fontsize = LabelSize)
plt.ylabel('System Response /nm (Arb units)', fontsize = LabelSize)
plt.xlabel("Wavelength (nm)", fontsize = LabelSize)
plt.legend(fontsize = LegSize)
plt.xlim(280,700)

plt.figure(10, figsize = (24,10))
plt.minorticks_on()
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
plt.scatter(CLData[0][0],CLData[0][1], linewidth=1, label = 'Raw Data', color = 'darkgrey')
plt.plot(CLData[0][0],yFilt10, linewidth=3, label = f'Savitzky-Golay Filter window = {window}, poly order = {poly},  R^2 = {np.around(r2_score(CLData[0][1],yFilt10),5)}', color = 'k', linestyle = '--')
plt.title('System response suggestion CL', fontsize = LabelSize)
plt.ylabel('System Response /nm (Arb units)', fontsize = LabelSize)
plt.xlabel("Wavelength (nm)", fontsize = LabelSize)
plt.legend(fontsize = LegSize)
plt.xlim(280,700)

plt.figure(11, figsize = (24,10))
plt.minorticks_on()
plt.xticks(fontsize=TickSize)
plt.yticks(fontsize=TickSize)
plt.scatter(CLData[0][0],CLData[0][1], linewidth=1, label = 'Raw Data', color = 'darkgrey')
plt.plot(CLData[0][0],yFilt11, linewidth=3, label = f'Whittaker-Eilers Filter Roughness penalty = {lam},\nSmoothing Order = {smoothOrder},  R^2 = {np.around(r2_score(CLData[0][1],yFilt11),5)}', color = 'k', linestyle = '--')
plt.title('System response suggestion CL', fontsize = LabelSize)
plt.ylabel('System Response /nm (Arb units)', fontsize = LabelSize)
plt.xlabel("Wavelength (nm)", fontsize = LabelSize)
plt.legend(fontsize = LegSize)
plt.xlim(280,700)

# np.savetxt('Output\\CLsystemResponse.txt',np.c_[CLData[0][0],yFilt11], delimiter=",", fmt='%.7f')