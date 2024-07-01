import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as matcol
import os
import csv
from scipy import stats
from scipy.optimize import curve_fit
import math
import Functions as F
from sklearn.metrics import r2_score
import nPlotterFunctions as nPF
import random

def CTI(value):
    return value/2.54

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

    plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))

DeltaLam = 0.3
h = 6.62607015e-34
c = 299792458

path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath('PLplotForIWGO.py')) + '\\PLTXT\\IWGO'))

for i in range(len(dirs)):
    vars()[f'files {i}'], vars()[f'data {i}'] = nPF.txtGrabber('PLplotForIWGO.py','\\PLTXT\\IWGO\\'+ dirs[i],True)

for i in range(len(dirs)):
    for j in range(len(vars()[f'files {i}'])):
        
        vars()[f'data {i}'][j][1] = vars()[f'data {i}'][j][1] * (((vars()[f'data {i}'][j][0]*1e-9)**2)/(h*c))
        vars()[f'data {i}'][j][0] = ((h*c) / (vars()[f'data {i}'][j][0]*1e-9)) * 6.242e18
        


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple','red','black','gold','orange','salmon','silver','pink','teal','magenta','darkgrey']
LineStyle = ['solid','dashed']

num = len(dirs)
plt.figure(1)
SubRow,SubCol = choose_subplot_dimensions(num)
plt.subplots(figsize=(CTI(16),CTI(16)),dpi = 600)
plt.suptitle('Photoluminescence of Ga₂O₃',fontsize = 14)

for i in range(len(dirs)):

    # add a new subplot iteratively
    plt.minorticks_on()
    plt.tight_layout()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(hspace=0.4,wspace=0)
    plt.xlim(1.75,4.5)
    ax = plt.subplot(SubRow, SubCol, i + 1)

    # filter df and plot ticker on the new subplot axis
    #RanCol = random.choice(list(matcol.XKCD_COLORS.values()))
    for j in range(len(vars()[f'files {i}'])):
        ax.plot(vars()[f'data {i}'][j][0], vars()[f'data {i}'][j][1], linestyle = LineStyle[j],color = LineColour[i], label = vars()[f'files {i}'][j])
    if i == int(np.round(num / 2)) - 1:
        ax.set_ylabel('Intensity (A.U)',fontsize = 14)
    else:
        pass
    
    ax.legend(loc = 'best')
    ax.set_title(LabList[i], fontsize = 14)

plt.minorticks_on()
plt.tight_layout()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.xlim(1.75,4.5)
plt.xlabel('Energy (eV)',fontsize = 14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


# LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']

# plt.figure(40,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# plt.minorticks_on()
# plt.grid(axis='both',which = 'major',linestyle=':', color = 'black')
# plt.grid(axis='both',which = 'minor',linestyle=':')
# plt.title('n values for different Ga₂O₃ polymorphs', fontsize = 12)
# plt.ylabel('n', fontsize = 12)
# plt.xlabel(r'$\lambda$ (nm)', fontsize = 12)
# for i in range(len(files)):
#         plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
#         plt.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
# for i in range(len(files1)):
#         plt.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
#         plt.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
# plt.legend()
# for i in range(len(files2)):
#         plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+2] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
#         plt.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+2],label = r'$n_{\kappa}$' + f'= {np.around(vars()[f"n2CoefList {files2[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files2[i]}"][0],3):.2e}/λ^2')
# plt.legend(fontsize = 16)