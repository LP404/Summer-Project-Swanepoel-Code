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
import matplotlib as mpl





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



h = 6.62607015e-34
c = 299792458


DeltaLam = 0.3

path, dirs, files, data, header = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\α','.csv')
path1, dirs1, files1, data1, header1 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\β','.csv')
path2, dirs2, files2, data2, header2 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\κ','.csv')

path3, dirs3, files3 = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\Absorption'))


fileNamesA, dataA = nPF.txtGrabber('BandGapCalc.py', '\\Data\\IWGO', '.txt', True)
thickness = np.array([350e-7,750e-7,740e-7])


for i in range(len(fileNamesA)):
    
    dataA[i][1] = dataA[i][1] /100
    
    vars()[f'Absorb 1{i}'] = np.log((1/dataA[i][1]))
    vars()[f'Absorb Fin{i}'] = vars()[f'Absorb 1{i}'] / thickness[i]
    
    dataA[i][0] = ((h*c) / (dataA[i][0]*1e-9)) * 6.242e18
    
    
    plt.plot(dataA[i][0],vars()[f'Absorb Fin{i}']**2 / max(vars()[f'Absorb Fin{i}']**2))
    plt.xlim(4.25,5.75)
    plt.ylim(0,0.6)
    

suffix = ".csv"

for i in range(len(files3)):
    files3[i] = files3[i][:-len(suffix)]  
    
    rows1 = []

    with open(str(path3)+"\\"+files3[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header1 = next(csvreader)
    
        for row1 in csvreader:
            rows1.append(row1)
    
    vars()[files3[i]+'_Data'] = rows1
    
    for j in range(len(vars()[files3[i]+'_Data'])):
        for k in range(len(vars()[files3[i]+'_Data'][j])):
                vars()[files3[i]+'_Data'][j][k] = float(vars()[files3[i]+'_Data'][j][k])


#Ind 0,5,6 n1
#Ind 0,13,14


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple','red','black','gold','orange','salmon','silver','pink','teal','magenta','darkgrey']

xP = np.linspace(190,800, 10001)

if data != 0:
    for i in range(len(files)):
        vars()[f"n1yFit {files[i]}"], vars()[f"n2yFit {files[i]}"], vars()[f"n1Cauch {files[i]}"], vars()[f"n2Cauch {files[i]}"], vars()[f"n1CoefList {files[i]}"], vars()[f"n2CoefList {files[i]}"] = nPF.CauchyFinder(files[i],data[i],xP,5,13,0)
    
    for i in range(len(files)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files[i]}"][k],vars()[f"n2yFit {files[i]}"][k]]
                
                writer.writerow(dataWrite)
else:
    pass

if data1 != 0:
    for i in range(len(files1)):
        vars()[f"n1yFit {files1[i]}"], vars()[f"n2yFit {files1[i]}"], vars()[f"n1Cauch {files1[i]}"], vars()[f"n2Cauch {files1[i]}"], vars()[f"n1CoefList {files1[i]}"], vars()[f"n2CoefList {files1[i]}"] = nPF.CauchyFinder(files1[i],data1[i],xP,5,13,0)

    for i in range(len(files1)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files1[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files1[i]}"][k],vars()[f"n2yFit {files1[i]}"][k]]
                
                writer.writerow(dataWrite)
else:
    pass

if data2 != 0:
    for i in range(len(files2)):
        vars()[f"n1yFit {files2[i]}"], vars()[f"n2yFit {files2[i]}"], vars()[f"n1Cauch {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n1CoefList {files2[i]}"], vars()[f"n2CoefList {files2[i]}"] = nPF.CauchyFinder(files2[i],data2[i],xP,5,13,0)
   
    for i in range(len(files2)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files2[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files2[i]}"][k],vars()[f"n2yFit {files2[i]}"][k]]
                
                writer.writerow(dataWrite)
else:
    pass

LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']

plt.figure(40,figsize=(CTI(16),CTI(16)), dpi=600)
plt.minorticks_on()
plt.grid(axis='both',which = 'major',linestyle=':', color = 'black')
plt.grid(axis='both',which = 'minor',linestyle=':')
plt.title('n values for different Ga₂O₃ polymorphs', fontsize = 14)
plt.ylabel('n', fontsize = 14)
plt.xlabel(r'$\lambda$ (nm)', fontsize = 14)

label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size 

label_size = 14
mpl.rcParams['ytick.labelsize'] = label_size 


for i in range(len(files)):
        plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        plt.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
for i in range(len(files1)):
        plt.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        plt.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
plt.legend()
for i in range(len(files2)):
        plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+2] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        plt.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+2],label = r'$n_{\kappa}$' + f'= {np.around(vars()[f"n2CoefList {files2[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files2[i]}"][0],3):.2e}/λ^2')
plt.legend(fontsize = 16)