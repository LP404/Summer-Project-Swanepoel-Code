import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

#Note the transmittance here is an idealised one with no frindges

def NoiseFilter(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    zi = lfilter_zi(b, a)
    z, _ = lfilter(b, a, yArray, zi=zi*yArray[0])
    z2, _ = lfilter(b, a, z, zi=zi*z[0])
    yFiltered = filtfilt(b, a, yArray)
    
    return yFiltered


def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]


path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\Absorption'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\Thickness'))


suffix = ".csv"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
    rows = []

    with open(str(path)+"\\"+files[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header = next(csvreader)
    
        for row in csvreader:
            rows.append(row)
    
    vars()[files[i]+'_Data'] = rows
    
    for j in range(len(vars()[files[i]+'_Data'])):
        for k in range(len(vars()[files[i]+'_Data'][j])):
                vars()[files[i]+'_Data'][j][k] = float(vars()[files[i]+'_Data'][j][k])

for i in range(len(files1)):
    files1[i] = files1[i][:-len(suffix)]  
    
    rows = []

    with open(str(path1)+"\\"+files1[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header1 = next(csvreader)
    
        for row in csvreader:
            rows.append(row)
    
    vars()[files1[i]+'_Data'] = rows
    
    for j in range(len(vars()[files1[i]+'_Data'])):
        for k in range(len(vars()[files1[i]+'_Data'][j])):
                vars()[files1[i]+'_Data'][j][k] = float(vars()[files1[i]+'_Data'][j][k])                
                
for i in range(len(files)):
    vars()[files[i]+'lambda'] = np.array(ListExtract(vars()[files[i]+'_Data'],0))
    vars()[files[i]+'Trans'] = np.array(ListExtract(vars()[files[i]+'_Data'],1))
    vars()[files[i]+'hv'] = np.array(ListExtract(vars()[files[i]+'_Data'],5))
    vars()[files[i]+'Alpha'] = np.array(ListExtract(vars()[files[i]+'_Data'],7))
    vars()[files[i]+'Alpha2'] = np.array(ListExtract(vars()[files[i]+'_Data'],9))
 
for i in range(len(files1)):
    vars()[files1[i]+'Thickness'] = np.array(ListExtract(vars()[files1[i]+'_Data'],2)) 

StepSize = 0.1

for i in range(len(files)):
    vars()[files[i]+'Trans'] = vars()[files[i]+'Trans'] * 100
    vars()[files[i]+'Trans'] = NoiseFilter(3,0.2,vars()[files[i]+'Trans'])
    vars()[files[i]+'DirectAllowed'] = (vars()[files[i]+'Alpha2']*vars()[files[i]+'hv'])**2
    vars()[files[i]+'IndirectAllowed'] = (vars()[files[i]+'Alpha2']*vars()[files[i]+'hv'])**0.5
    vars()[files[i]+'DirectDisallowed'] = (vars()[files[i]+'Alpha2']*vars()[files[i]+'hv'])**(2/3)
    vars()[files[i]+'IndirectDisallowed'] = (vars()[files[i]+'Alpha2']*vars()[files[i]+'hv'])**(1/3)
    
    
    vars()[files[i]+'ticks'] = np.arange(min(vars()[files[i]+'hv']),max(vars()[files[i]+'hv'])+StepSize,StepSize)

for i in range(len(files)):
    vars()[files[i]+'NewAbsorb'] = ((np.log((1/vars()[files[i]+'Trans']))) / (vars()[files1[i]+'Thickness'] * 1e-9))

for i in range(len(files)):
    plt.figure(i,figsize=(15,12))
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'lambda'],vars()[files[i]+'Trans'])
    plt.title(files[i]+'_Transmittance')
    plt.ylabel('Transmittance')
    plt.xlabel('hv (eV)')
    # plt.xticks(vars()[files[i]+'ticks'],rotation = 90,fontsize=5)
    plt.plot()
    
    plt.figure(i+len(files)+1,figsize=(15,12))
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.title(files[i]+'_Absorbance')
    plt.plot(vars()[files[i]+'hv'],vars()[files[i]+'Alpha2'])
    plt.ylabel('α (cm^-1)')
    plt.xlabel('hv (eV)')
    # plt.xticks(vars()[files[i]+'ticks'],rotation = 90,fontsize=5)
    plt.plot()
    
    plt.figure(i+len(files)+1*(len(files))+1,figsize=(15,12), dpi =1200)
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.title(files[i]+'Direct Allowed')
    plt.plot(vars()[files[i]+'hv'],vars()[files[i]+'DirectAllowed'], label = 'Direct Allowed')
    #plt.plot(vars()[files[i]+'hv'],vars()[files[i]+'IndirectAllowed'], label = 'Indirect Allowed')
    plt.ylabel('α (cm^-1)')
    plt.xlabel('hv (eV)')
    plt.legend()
    plt.xticks(vars()[files[i]+'ticks'],rotation = 90,fontsize=5)
    plt.plot()
    


