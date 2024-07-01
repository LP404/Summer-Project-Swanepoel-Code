import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import csv
from scipy import stats
from scipy.optimize import curve_fit
import math
import Functions as F
from sklearn.metrics import r2_score
import nPlotterFunctions as nPF
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from decimal import Decimal
import random
import matplotlib as mpl
import nPlotterFunctions as nPF
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi
from natsort import natsorted
from scipy.ndimage import gaussian_filter1d

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


def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]


def choose_subplot_dimensions(k):
    if k < 4:
        return k, 1
    elif k < 11:
        return math.ceil(k/2), 2
    else:
        # I've chosen to have a maximum of 3 columns
        return math.ceil(k/3), 3

def TickFinder(val):
    if val >= 10:
        return int(np.log10(val))
    else:
        var = float(Decimal(val) % 1)
        if str(val)[::-1].find('.') < str(var)[::-1].find('.'):
            return np.around(var,(str(val)[::-1].find('.')+1))
        else:
            return var

def ExtractIntFromStr(Str):
    return [int(s) for s in Str.split() if s.isdigit()][0]

def underline_annotation(text):
    f = plt.gcf()
    ax = plt.gca()
    tb = text.get_tightbbox(f.canvas.get_renderer()).transformed(f.transFigure.inverted())
                            # text isn't drawn immediately and must be 
                            # given a renderer if one isn't cached.
                                                    # tightbbox return units are in 
                                                    # 'figure pixels', transformed 
                                                    # to 'figure fraction'.
    ax.annotate('', xy=(tb.xmin,tb.y0), xytext=(tb.xmax,tb.y0),
                xycoords="figure fraction",
                arrowprops=dict(arrowstyle="-", color='k'))
    
    plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))


axesLabelSize = 32
axesLabelPad = 32
axesFontSize = 26 
axesFontPadding = 20
lineWidth = 5
legendFontSize = 24

GraphResolution = 100
sizeFigSqaure = CTI(30)

DeltaLam = 0.3

tickMajorSize = 15
tickMinorSize = 10


plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['blue', 'forestgreen', 'darkred', 'cyan', 'magenta', 'black', 'purple', 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold'])

# path, dirs, files, data = nPF.fileGrabberCSVNH('GenericCSVPlotter.py','\\BiOI\\Input\\03_2024\\CSV1\\Left')
# path, dirs, files, data = nPF.fileGrabberCSVNH('GenericCSVPlotter.py','\\BiOI\\Input\\03_2024\\CSV1\\Middle')
# path, dirs, files, data = nPF.fileGrabberCSVNH('GenericCSVPlotter.py','\\BiOI\\Input\\03_2024\\CSV1\\Right')

# path, dirs, files, data = nPF.fileGrabberCSVNH('GenericCSVPlotter.py','\\BiOI\\Input\\03_2024\\CSV2\\Left')
path, dirs, files, data = nPF.fileGrabberCSVNH('GenericCSVPlotter.py','\\BiOI\\Input\\03_2024\\CSV2\\Right')



for i in range(len(files)):
    vars()[f'x {[i]}'] = np.array(ListExtract(data[i],0))
    vars()[f'y {[i]}'] = np.array(ListExtract(data[i],1))


# for i in range(len(files)):
#     vars()[f'x {[i]}'], vars()[f'y {[i]}'] = ConvertDependance(np.array(ListExtract(data[i],0)),np.array(ListExtract(data[i],1)))

#     np.flip(vars()[f'x {[i]}'])
#     np.flip(vars()[f'y {[i]}'])
    
for i in range(len(files)):
    plt.figure(i,figsize=(CTI(20),CTI(10)), dpi=600)
    plt.minorticks_on()
    plt.grid(axis='both',which = 'major',linestyle=':', color = 'black')
    plt.grid(axis='both',which = 'minor',linestyle=':')
    plt.title(f'{files[i]} spectrum at 80K', fontsize = 12)
    plt.ylabel('Intensity (Arb. Units)', fontsize = 12)
    plt.xlabel('Energy (eV)', fontsize = 12)
    plt.plot(vars()[f'x {[i]}'], gaussian_filter1d(vars()[f'y {[i]}'],20), linewidth = 1)
    plt.xlim(min(vars()[f'x {[i]}']),2.2)
    
HowManyNotPoint = 0
GSTD = 20

plt.figure(len(files)+1,figsize=(CTI(20),CTI(10)), dpi=600)
plt.minorticks_on()
plt.grid(axis='both',which = 'major',linestyle=':', color = 'black')
plt.grid(axis='both',which = 'minor',linestyle=':')
plt.title('Region Intensity spectrum at 80K', fontsize = 12)
plt.ylabel('Intensity (Arb. Units)', fontsize = 12)
plt.xlabel('Energy (eV)', fontsize = 12)
for i in range(HowManyNotPoint,len(files)):
# for i in range(HowManyNotPoint,4):
     plt.plot(vars()[f'x {[i]}'], vars()[f'y {[i]}'], linewidth = 1, label = f'{files[i]}')
#    plt.plot(vars()[f'x {[i]}'], (gaussian_filter1d(vars()[f'y {[i]}'],GSTD) / max(gaussian_filter1d(vars()[f'y {[i]}'],GSTD))), linewidth = 1, label = f'{files[i]}')
plt.xlim(min(vars()[f'x {[i]}']),2.2)
handles, labels = plt.gca().get_legend_handles_labels()
order = natsorted(range(len(files)), key=lambda i: files[i])
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],fontsize = 8)
#plt.legend(fontsize = 8)