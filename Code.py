import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import Functions as F

#dt = 1


#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\ThorLabsData'))

xMin = 200
xMax = 875

#Provisional will change later
nSub = 1.75

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i].rstrip(".txt")
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)

    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)], vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],10)

    vars()['MaxBlocks'+str(i)] = F.BlockFinder(vars()['MaxX'+str(i)])
    vars()['MinBlocks'+str(i)] = F.BlockFinder(vars()['MinX'+str(i)])

    vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)] = F.AntiNodeHighlander(vars()['MaxBlocks'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])
    vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)] = F.AntiNodeHighlander(vars()['MinBlocks'+str(i)],vars()['MinX'+str(i)], vars()['MinY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])

    #Put case statment here in python 3.10
    if i == 0:
        #All the xValues are the same, only need to do this once
        xP = np.linspace(min(vars()[files[i]+'T'][0]),max(vars()[files[i]+'T'][0]),10001)
    else:
        pass

    vars()['yPMax'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)])
    vars()['yPMin'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)])


for i in range(len(files1)):
    files1[i] = files1[i].rstrip(".txt")
    
    vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files1[i]][0] = vars()[files1[i]][0] * 1000
    vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
    vars()[files1[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files1[i]],xMin,xMax))

    vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])

for i in range(len(files)):
    plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1])
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    plt.scatter(vars()['xNewMax'+str(i)] , vars()['yNewMax'+str(i)], color = 'black', marker = "x")
    plt.scatter(vars()['xNewMin'+str(i)] , vars()['yNewMin'+str(i)], color = 'red', marker = "x")
    plt.plot(xP,vars()['yPMax'+str(i)], color = 'black', linestyle="dotted")
    plt.plot(xP,vars()['yPMin'+str(i)], color = 'red', linestyle="dotted")
    
    
