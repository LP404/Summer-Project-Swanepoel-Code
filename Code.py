import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import Functions as F

#dt = 1


#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
xMin = 200
xMax = 875

#Provisional will change later
nSub = 1.75

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i].rstrip(".txt")

#The next three loops import and assign varable names to the data 
for i in range(len(files)):
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T

for i in range(len(files)):
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]][0],vars()[files[i]][1],xMin,xMax)

#T stands for Truncated
    
    
for i in range(len(files)):
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)], vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],10)

for i in range(len(files)):
    vars()['MaxBlocks'+str(i)] = F.BlockFinder(vars()['MaxX'+str(i)])
    vars()['MinBlocks'+str(i)] = F.BlockFinder(vars()['MinX'+str(i)])

for i in range(len(files)):
    vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)] = F.AntiNodeHighlander(vars()['MaxBlocks'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])
    vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)] = F.AntiNodeHighlander(vars()['MinBlocks'+str(i)],vars()['MinX'+str(i)], vars()['MinY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])

#All the xValues are the same, only need to do this once
xP = np.linspace(min(vars()[files[0]+'T'][0]),max(vars()[files[0]+'T'][0]),10001)

for i in range(len(files)):
    vars()['yPMax'+str(i)], vars()['MaxLocWave'+str(i)] = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)])
    vars()['yPMin'+str(i)], vars()['MinLocWave'+str(i)] = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)])






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
    
    
