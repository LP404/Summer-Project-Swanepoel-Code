import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import copy
import random
from whittaker_eilers import WhittakerSmoother
from sklearn.metrics import r2_score
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi, savgol_filter
from scipy.ndimage import gaussian_filter1d
from sklearn.metrics import r2_score
import pandas as pd
from natsort import natsorted
import matplotlib.patches as mpatches

def xlsxImport(filePath):
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(__file__)) + filePath))

    fileNames = natsorted(fileNames)
    
    data = []
    
    
    for i in range(len(fileNames)):
        fileNames[i] = fileNames[i][:-5]
        
        vars()['df_'+fileNames[i]] = pd.read_excel(path + "\\" + fileNames[i] + ".xlsx", header=None).T.to_numpy()
        
        data.append(vars()['df_'+fileNames[i]])

    
    
    return fileNames, data

def IntInt(x,y,genCase):
    
    totalX = 0
    
    for i in range(1,len(x)):
        
        totalX += abs(x[i]-x[i-1])
    
    meanDelX = totalX / len(x)
    
    if genCase == True:
        IntIntensity = sum(y) * meanDelX    
    else:
        IntIntensity = sum(y) * 0.0007717666774969188
    
    return IntIntensity

def nListMaker(n):
    return [[] for x in range(n)]


def BinNoSelector(maxVal,width):
    
    BinNo = maxVal/width
    
    return int(BinNo)
    
    
    


RefVar = ['3CA1','3EA1','3CA2','3EA2','TCA1','TEA1','TCA2','TEA2']

ISC = [True,False,True,False,True,False,True,False]
IS3 = [True,True,True,True,False,False,False,False]
ISA1 = [True,True,False,False,True,True,False,False]


fPath = ['\\3\\Centre\\Area1',
         '\\3\\Edge\\Area1',
         '\\3\\Centre\\Area2',
         '\\3\\Edge\\Area2',
         '\\Treated\\Centre\\Area1',
         '\\Treated\\Edge\\Area1',
         '\\Treated\\Centre\\Area2',
         '\\Treated\\Edge\\Area2']


title = ['Sample 3 Centre Spectra: Area 1',
         'Sample 3 Edge Spectra: Area 1',
         'Sample 3 Centre Spectra: Area 2',
         'Sample 3 Edge Spectra: Area 2',
         'Treated Centre Spectra: Area 1',
         'Treated Edge Spectra: Area 1',
         'Treated Centre Spectra: Area 2',
         'Treated Edge Spectra: Area 2']


title1 = ['Sample 3 Centre Integrated Intensity: Area 1',
          'Sample 3 Edge Integrated Intensity: Area 1',
          'Sample 3 Centre Integrated Intensity: Area 2',
          'Sample 3 Edge Integrated Intensity: Area 2',
          'Treated Centre Integrated Intensity: Area 1',
          'Treated Edge Integrated Intensity: Area 1',
          'Treated Centre Integrated Intensity: Area 2',
          'Treated Edge Integrated Intensity: Area 2']
  

for i in range(len(RefVar)):
    
    fName, data  = xlsxImport(fPath[i])
    
    vars()[RefVar[i] + '_data'] = data
    
    lam = 3e5
    smoothOrder = 3
    
    vars()['intint_' + RefVar[i]] = np.array([])
    
    
    plt.figure(i)
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    
    plt.title(title[i])
    plt.ylabel('Intensity (counts/eV)')
    plt.xlabel('Energy  (eV)')

    for j in range(len(vars()[RefVar[i] + '_data'])):
        
        ws = WhittakerSmoother(lmbda=lam, order=smoothOrder, data_length = len(vars()[RefVar[i] + '_data'][j][1]))
        y1 = ws.smooth(vars()[RefVar[i] + '_data'][j][1])
        y2= ws.smooth_optimal(vars()[RefVar[i] + '_data'][j][1], break_serial_correlation=False)
        optimally_smoothed_series = y2.get_optimal().get_smoothed()
        optimal_lambda = y2.get_optimal().get_lambda()
        
        if ISA1[i] == True: 
            plt.plot(vars()[RefVar[i] + '_data'][j][0], np.array(optimally_smoothed_series) / max(optimally_smoothed_series), label = fName[j])
        else:
            plt.plot(vars()[RefVar[i] + '_data'][j][0], np.array(optimally_smoothed_series) / max (optimally_smoothed_series), label = str(int(fName[j]) - 10))
        
        vars()['intint_' + RefVar[i]] = np.append(vars()['intint_' + RefVar[i]],IntInt(vars()[RefVar[i] + '_data'][j][0],vars()[RefVar[i] + '_data'][j][1],False))
        
        print(np.where(vars()[RefVar[i] + '_data'][j][1] == max(vars()[RefVar[i] + '_data'][j][1]))[0][0])
        
    plt.legend()


### ISC
### IS3
### ISA1


IntbinCent3 = np.array([]) 
IntbinEdge3 = np.array([]) 
    
IntbinCentT = np.array([]) 
IntbinEdgeT = np.array([]) 


for i in range(len(RefVar)):
    
    plt.figure(100)
    plt.title('Center And Grain spectra for treated sample')
    plt.ylabel('Intensity (counts/eV)')
    plt.xlabel('Energy  (eV)')
    plt.ylim(0,3.65e6)
    plt.xlim(1.25482782109334,2.2)
    if ISC[i] == True and IS3[i] == False:
        for j in range(len(data)):
            plt.plot(vars()[RefVar[i] + '_data'][j][0], vars()[RefVar[i] + '_data'][j][1], color ='grey')
    elif ISC[i] == False and IS3[i] == False:
        for j in range(len(data)):
            plt.plot(vars()[RefVar[i] + '_data'][j][0], vars()[RefVar[i] + '_data'][j][1],color= 'red', alpha=0.4) 

    red_patch = mpatches.Patch(color='red', label='Edge', alpha = 0.4)
    grey_patch = mpatches.Patch(color='grey', label='Centre')
     
    plt.legend(handles=[red_patch, grey_patch])
            
    
    plt.figure(101)
    plt.title('Center And Grain spectra for sample 3')
    plt.ylabel('Intensity (counts/eV)')
    plt.xlabel('Energy  (eV)')
    plt.ylim(0,2.5e6)
    plt.xlim(1.25482782109334,2.2)
    if ISC[i] == True and IS3[i] == True:
        for j in range(len(data)):
            plt.plot(vars()[RefVar[i] + '_data'][j][0], vars()[RefVar[i] + '_data'][j][1], color ='grey')
    elif ISC[i] == False and IS3[i] == True:
        for j in range(len(data)):
            plt.plot(vars()[RefVar[i] + '_data'][j][0], vars()[RefVar[i] + '_data'][j][1],color= 'red', alpha=0.4)   

    red_patch = mpatches.Patch(color='red', label='Edge', alpha = 0.4)
    grey_patch = mpatches.Patch(color='grey', label='Centre')
     
    plt.legend(handles=[red_patch, grey_patch])

    # plt.figure(102)
    # plt.title('Center And Grain integrated intensites for treated sample')
    # plt.ylabel('Counts ')
    # plt.xlabel('Integrated intensity')
    # plt.ylim(0,3.65e6)
    # plt.xlim(1.25482782109334,2.2)

    
    # np.concatenate((arr5,arr4),axis=None)
    
    ##Check bin lengths
    plt.figure(102)
    plt.title('Histogram of Integrated Intensites of Sample 3')
    plt.xlabel('Integrated Intensity (Arb.)')
    plt.ylabel('Occurrences')
    if ISC[i] == True and IS3[i] == True:
        IntbinCent3 = np.concatenate((IntbinCent3,vars()['intint_' + RefVar[i]]),axis=None)
        
        
        if ISA1[i] == True:
            pass        
        else:
            plt.hist(IntbinCent3, color = 'grey', label = 'Centre', bins = BinNoSelector(max(IntbinCent3),50000))
    
    elif ISC[i] == False and IS3[i] == True:
        IntbinEdge3 = np.concatenate((IntbinEdge3,vars()['intint_' + RefVar[i]]),axis=None)
        
        if ISA1[i] == True:  
            pass
        else:
            plt.hist(IntbinEdge3, alpha= 0.5, color = 'red', label = 'Edge', bins = BinNoSelector(max(IntbinEdge3),50000))    
    
    red_patch = mpatches.Patch(color='red', label= f'Edge: μ = {np.around(np.mean(IntbinEdge3),2)} ± {np.around((np.std(IntbinEdge3) / (20**0.5)),2)} ', alpha = 0.5)
    grey_patch = mpatches.Patch(color='grey', label= f'Centre: μ = {np.around(np.mean(IntbinCent3),2)} ± {np.around((np.std(IntbinCent3) / (20**0.5)),2)}')
    white_patch = mpatches.Patch(color='white', label='Bin width ≈ 50,000 units', alpha = 0)
     
    
    plt.legend(handles=[red_patch, grey_patch, white_patch])

    plt.figure(103)
    plt.title('Histogram of Integrated Intensites of Treated Sample')
    plt.xlabel('Integrated Intensity (Arb.)')
    plt.ylabel('Occurrences')
    if ISC[i] == True and IS3[i] == False:
        IntbinCentT = np.concatenate((IntbinCentT,vars()['intint_' + RefVar[i]]),axis=None)

        
        if ISA1[i] == True:
            pass
        else:             
            plt.hist(IntbinCentT, color = 'grey', label = 'Centre', bins = BinNoSelector(max(IntbinCentT),50000))
            
    elif ISC[i] == False and IS3[i] == False:
        IntbinEdgeT = np.concatenate((IntbinEdgeT,vars()['intint_' + RefVar[i]]),axis=None)

        
        if ISA1[i] == True:
            pass
        else:
            plt.hist(IntbinEdgeT, alpha= 0.5, color = 'red', label = 'Edge', bins = BinNoSelector(max(IntbinEdgeT),50000))
    
    red_patch = mpatches.Patch(color='red', label= f'Edge: μ = {np.around(np.mean(IntbinEdgeT),2)} ± {np.around((np.std(IntbinEdgeT) / (20**0.5)),2)} ', alpha = 0.5)
    grey_patch = mpatches.Patch(color='grey', label= f'Centre: μ = {np.around(np.mean(IntbinCentT),2)} ± {np.around((np.std(IntbinCentT) / (20**0.5)),2)} ')
    white_patch = mpatches.Patch(color='white', label='Bin width ≈ 50,000 units', alpha = 0)
     
    plt.legend(handles=[red_patch, grey_patch, white_patch])
        
   
    
    # plt.title('Histogram of Integrated Intensites')
    # plt.xlabel('Int intensites')
    # plt.ylabel('Counts')
    # plt.hist(arr1, color = 'grey', label = 'Centre', bins=16)
    # plt.hist(arr2, color = 'darkred', alpha = 0.5,label = 'Edge', bins=10)
    # plt.legend()
