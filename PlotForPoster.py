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
from decimal import Decimal

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


def sellmerier(lam,b1,b2,b3,c1,c2,c3):
    return np.sqrt(( 1 + ((b1*lam**2)/(lam**2 - c1)) + ((b2*lam**2)/(lam**2 - c2)) + ((b3*lam**2)/(lam**2 - c3))))


def TickFinder(val):
    if val >= 10:
        return int(np.log10(val))
    else:
        var = float(Decimal(val) % 1)
        if str(val)[::-1].find('.') < str(var)[::-1].find('.'):
            return np.around(var,(str(val)[::-1].find('.')+1))
        else:
            return var


    
plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))




DeltaLam = 0.3

path, dirs, files, data, header = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\α','.csv')
path1, dirs1, files1, data1, header1 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\β','.csv')
path2, dirs2, files2, data2, header2 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\κ','.csv')


pathA, dirsA, filesA, dataA, headerA = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\α','.csv')
path1A, dirs1A, files1A, data1A, header1A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\β','.csv')
path2A, dirsA, files2A, data2A, header2A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\κ','.csv')

pathT, dirsT, filesT, dataT, headerT = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\α','.csv')
path1T, dirs1T, files1T, data1T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\β','.csv')
path2T, dirs2T, files2T, data2T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\κ','.csv')

suffix = ".csv"

#Ind 0,5,6 n1
#Ind 0,13,14


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','red','black','gold','orange','salmon','silver','pink','teal','purple']

xP = np.linspace(190,800, 10001)

if data != 0:
    for i in range(len(files)):
        vars()[f"n1yFit {files[i]}"], vars()[f"n2yFit {files[i]}"], vars()[f"n1Cauch {files[i]}"], vars()[f"n2Cauch {files[i]}"], vars()[f"n1CoefList {files[i]}"], vars()[f"n2CoefList {files[i]}"] = nPF.CauchyFinder(files[i],data[i],xP,5,13,0)   
else:
    pass

if data1 != 0:
    for i in range(len(files1)):
        vars()[f"n1yFit {files1[i]}"], vars()[f"n2yFit {files1[i]}"], vars()[f"n1Cauch {files1[i]}"], vars()[f"n2Cauch {files1[i]}"], vars()[f"n1CoefList {files1[i]}"], vars()[f"n2CoefList {files1[i]}"] = nPF.CauchyFinder(files1[i],data1[i],xP,5,13,0)
else:
    pass

if data2 != 0:
    for i in range(len(files2)):
        vars()[f"n1yFit {files2[i]}"], vars()[f"n2yFit {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2CoefList {files2[i]}"], vars()[f"n2CoefList {files2[i]}"] = nPF.CauchyFinder(files2[i],data2[i],xP,5,13,0)

else:
    pass


fig = plt.figure(40, figsize = (4,4),dpi = 1200)
ax = fig.add_axes([0,0,1,1])
ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(xP))))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(xP)) / 4))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f"n2yFit {files[0]}"])),2) / 4))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f"n2yFit {files[0]}"])),2) / 16))
ax.set_title('n values of Ga₂O₃ polymorphs', pad=10)
ax.set_ylabel('n (arb. units)', labelpad=10)
ax.set_xlabel(r'$\lambda$ (nm)', labelpad=10)
for i in range(len(files)):
        ax.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        ax.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
for i in range(len(files1)):
        ax.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        ax.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
ax.legend(bbox_to_anchor=(0.95, 0.95), fontsize=10, title = 'Cauchy Equations')

