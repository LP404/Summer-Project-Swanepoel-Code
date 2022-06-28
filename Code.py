import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import Functions as F

#dt = 1


#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))

xMin = 200
xMax = 875

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

    vars()['nS'+files1[i]] = F.Sub_nFinder(vars()['yP'+files1[i]])
    
    

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
    
for i in range(len(files1)):   
    plt.plot(xP,vars()['yP'+files1[i]], color = 'orange')
    
   
for i in range(len(files)):
    
    vars()['xNewMaxT'+str(i)], vars()['yNewMaxT'+str(i)],vars()['xNewMinT'+str(i)], vars()['yNewMinT'+str(i)] = F.SecondTrim( vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)])
    
    vars()['dForMax'+str(i)] = np.array([])
    vars()['dForMin'+str(i)] = np.array([])
    vars()['nForMax'+str(i)] = np.array([])
    vars()['nForMin'+str(i)] = np.array([])
    vars()['mForMax'+str(i)] = np.array([])
    vars()['mForMin'+str(i)] = np.array([])
    
    for k in range(len(vars()['xNewMaxT'+str(i)])):
    
        
        Loc = np.where(xP == F.FindNearest(xP,vars()['xNewMaxT'+str(i)][k]))[0][0]
        Opp = vars()['yPMin'+str(i)][Loc]   
        Sub = vars()['nS'+files1[0]][Loc]
                
        n = F.nFinder(vars()['yNewMaxT'+str(i)][k],Opp,vars()['xNewMaxT'+str(i)][k],Sub)

        vars()['nForMax'+str(i)] = np.append(vars()['nForMax'+str(i)],n)
    
    for j in range(len(vars()['xNewMinT'+str(i)])):
        
        
        Loc1 = np.where(xP == F.FindNearest(xP,vars()['xNewMinT'+str(i)][j]))[0][0]
        Opp1 = vars()['yPMax'+str(i)][Loc1]   
        Sub1 = vars()['nS'+files1[0]][Loc1]
                
        n1 = F.nFinder(Opp1,vars()['yNewMinT'+str(i)][j],vars()['xNewMinT'+str(i)][j],Sub1)    
  
        vars()['nForMin'+str(i)] = np.append(vars()['nForMin'+str(i)],n1)
  
    
    for k in range(1,len(vars()['xNewMaxT'+str(i)])):
        
        d = F.dFinder(vars()['nForMax'+str(i)][k-1],vars()['nForMax'+str(i)][k],vars()['xNewMaxT'+str(i)][k-1],vars()['xNewMaxT'+str(i)][k])
    
        vars()['dForMax'+str(i)] = np.append(vars()['dForMax'+str(i)],abs(d))
    
    for j in range(1,len(vars()['xNewMinT'+str(i)])):
            
        d1 = F.dFinder(vars()['nForMin'+str(i)][j-1],vars()['nForMin'+str(i)][j],vars()['xNewMinT'+str(i)][j-1],vars()['xNewMinT'+str(i)][j])
    
        vars()['dForMin'+str(i)] = np.append(vars()['dForMin'+str(i)],abs(d1))
        
    
    for k in range(len(vars()['xNewMaxT'+str(i)])):
        if k == 0:
            mCalc =  2*vars()['nForMax'+str(i)][0]*vars()['dForMax'+str(i)][0] / vars()['xNewMaxT'+str(i)][0]
            vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
        elif k == (len(vars()['xNewMaxT'+str(i)])-1):
            mCalc =  2*vars()['nForMax'+str(i)][k]*vars()['dForMax'+str(i)][k-1] / vars()['xNewMaxT'+str(i)][k]
            vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
        else:
            for l in range(0,len(vars()['dForMax'+str(i)]) - 1):
                mCalc =  2*vars()['nForMax'+str(i)][k]*vars()['dForMax'+str(i)][k-1+l] / vars()['xNewMaxT'+str(i)][k]
                vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
    
    for j in range(len(vars()['xNewMinT'+str(i)])):
        if j == 0:
            mCalc1 =  2*vars()['nForMin'+str(i)][0]*vars()['dForMin'+str(i)][0] / vars()['xNewMinT'+str(i)][0]
            vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))
        elif j == (len(vars()['xNewMinT'+str(i)])-1):
            mCalc1 =  2*vars()['nForMin'+str(i)][j]*vars()['dForMin'+str(i)][j-1] / vars()['xNewMinT'+str(i)][j]
            vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))
        else:
            for o in range(0,len(vars()['dForMin'+str(i)]) - 1):
                mCalc1 =  2*vars()['nForMin'+str(i)][j]*vars()['dForMin'+str(i)][j-1+o] / vars()['xNewMinT'+str(i)][j]
                vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))



    # x1 = F.xFinder(n1,1.78,62.3486/100,59.6286/100)
    # x2 = F.xFinder(n2,1.775,63.6930/100,60.8116/100)

    # k1 = (-446e-9 / (4 * np.pi * d)) * np.log(x1)
    # k2 = (-534e-9 / (4 * np.pi * d)) * np.log(x2)

