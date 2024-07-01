import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import Functions1 as F1
import nPlotterFunctions as nPF
import csv
from scipy.optimize import curve_fit


SuperPath,SuperDirs,SuperFiles = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotterAgain.py')) + '\\DataComparison\\Samples'))

for i in range(len(SuperDirs)):
    vars()[f'path {i}'], vars()[f'dirs {i}'], vars()[f'files {i}'] = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotterAgain.py')) + '\\DataComparison\\Samples\\'+str(SuperDirs[i])))


pathA, dirsA, filesA = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotterAgain.py')) + '\\DataComparison\\Substrates\\Single'))
pathB, dirsB, filesB = next(os.walk(os.path.dirname(os.path.realpath('GenerticTxtPlotterAgain.py')) + '\\DataComparison\\Substrates\\Double'))

xMin = 200
xMax = 800

sizeFigSqaure = F1.CTI(20)

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
suffix = ".txt"
for i  in range(len(SuperDirs)):
    for j in range(len( vars()[f'files {i}'])):
        vars()[f'files {i}'][j] = vars()[f'files {i}'][j][:-len(suffix)]  
        
        vars()[f'files data {i}{j}'] = np.loadtxt(open(vars()[f'path {i}'] + "\\" + vars()[f'files {i}'][j] + ".txt", "rb"), delimiter=",",skiprows = 2).T
        vars()[f'files data {i}{j}'][1] = vars()[f'files data {i}{j}'][1] / 100
    
        #T stands for Truncated
        vars()[f'files data {i}{j} T'] = F.Trim(vars()[f'files data {i}{j}'],xMin,xMax)
    
    
        # #Change 0.075 back to 0.2
        # vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[f'files data {i}{j} T'][1])
    
    
    if i == 0:
        #All the xValues are the same, only need to do this once
        xP = np.linspace(min(vars()[f'files data {i}{j} T'][0]),max(vars()[f'files data {i}{j} T'][0]),10001)
    else:
        pass
    


for i in range(len(filesA)):
    filesA[i] = filesA[i].rstrip(".txt")
    
    vars()[filesA[i]] = np.loadtxt(open(pathA + "\\" + filesA[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    

    vars()[filesA[i]][0] = vars()[filesA[i]][0]
    
    vars()[filesA[i]][1] = vars()[filesA[i]][1] / 100
    

    vars()[filesA[i]+'T'] = F.Trim(vars()[filesA[i]],xMin,xMax)

    vars()['yP'+filesA[i]] = np.interp(xP,vars()[filesA[i]+'T'][0],vars()[filesA[i]+'T'][1])

for i in range(len(filesB)):
    filesB[i] = filesB[i].rstrip(".txt")
    
    vars()[filesB[i]] = np.loadtxt(open(pathB + "\\" + filesB[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    

    vars()[filesB[i]][0] = vars()[filesB[i]][0]
    
    vars()[filesB[i]][1] = vars()[filesB[i]][1] / 100
    

    vars()[filesB[i]+'T'] = F.Trim(vars()[filesB[i]],xMin,xMax)

    vars()['yP'+filesB[i]] = np.interp(xP,vars()[filesB[i]+'T'][0],vars()[filesB[i]+'T'][1])



# plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
# plt.figure(0,figsize = (8,5),dpi = 600)
for i  in range(len(SuperDirs)):
    plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.title(f'{SuperDirs[i]}')
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    for j in range(len( vars()[f'files {i}'])):

        plt.plot(vars()[f'files data {i}{j} T'][0],vars()[f'files data {i}{j} T'][1],label = vars()[f'files {i}'][j])
    # plt.plot(xP,vars()['yP'+filesA[1]], label = filesA[1])
    # plt.plot(xP,vars()['yP'+filesB[3]], label = filesB[3])
    plt.legend()


plt.figure(100,figsize=(29.7/2.54,21.0/2.54), dpi=600)
plt.minorticks_on()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.title('Sapphire Single Side Polish')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission/Absorbance")
for i in range(len(filesA)):
    plt.plot(xP,vars()['yP'+filesA[i]], label = filesA[i])
plt.legend()



plt.figure(200,figsize=(29.7/2.54,21.0/2.54), dpi=600)    
plt.minorticks_on()
plt.grid(which='major', color='k', linestyle='-')
plt.grid(which='minor', color='darkgray', linestyle='--')
plt.title('Sapphire Double Side Polish')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission/Absorbance")
for i in range(len(filesB)):
    plt.plot(xP,vars()['yP'+filesB[i]], label = filesB[i])
plt.legend()

