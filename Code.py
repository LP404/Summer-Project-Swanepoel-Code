import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import Functions as F
import pandas as pd

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

    vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)] = F.AntiNodeHighlander(vars()['MaxX'+str(i)], vars()['MaxY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])
    vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)] = F.AntiNodeHighlander(vars()['MinX'+str(i)], vars()['MinY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])

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
    plt.plot(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1], label = "Raw Data")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    
    #np.savetxt(str(files[i])+".csv", vars()[files[i]+'T'].T, delimiter=',')


    plt.scatter(vars()['xNewMax'+str(i)] , vars()['yNewMax'+str(i)], color = 'black', marker = "x", label = "Maxima")
    plt.scatter(vars()['xNewMin'+str(i)] , vars()['yNewMin'+str(i)], color = 'red', marker = "x", label = "Minima")
    plt.plot(xP,vars()['yPMax'+str(i)], color = 'black', linestyle="dotted", label = "TM")
    plt.plot(xP,vars()['yPMin'+str(i)], color = 'red', linestyle="dotted", label = "Tm")
    
for i in range(len(files1)):   
    plt.plot(xP,vars()['yP'+files1[i]], color = 'orange', label = "Substrate")
    
plt.legend()
    

for i in range(len(files)):
    
    vars()['xNewMaxT'+str(i)], vars()['yNewMaxT'+str(i)],vars()['xNewMinT'+str(i)], vars()['yNewMinT'+str(i)] = F.SecondTrim( vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)])
    
    vars()['TmForMax'+str(i)] = np.array([])
    vars()['TMForMin'+str(i)] = np.array([])    
    vars()['dForMax'+str(i)] = np.array([])
    vars()['dForMin'+str(i)] = np.array([])
    vars()['nForMax'+str(i)] = np.array([])
    vars()['nForMin'+str(i)] = np.array([])
    vars()['mForMax'+str(i)] = np.array([])
    vars()['mForMin'+str(i)] = np.array([])
    vars()['New_dForMax'+str(i)] = np.array([])
    vars()['New_dForMin'+str(i)] = np.array([])
    vars()['New_nForMax'+str(i)] = np.array([])
    vars()['New_nForMin'+str(i)] = np.array([])
    
    HasLoopMax = 0
    HasLoopMin = 0
    
    HasLoopMax1 = 0
    HasLoopMin1 = 0
    
    HasLoopMax2 = 0
    HasLoopMin2 = 0
    
    for k in range(len(vars()['xNewMaxT'+str(i)])):
    
        
        Loc = np.where(xP == F.FindNearest(xP,vars()['xNewMaxT'+str(i)][k]))[0][0]
        Opp = vars()['yPMin'+str(i)][Loc]   
        vars()['TmForMax'+str(i)] = np.append(vars()['TmForMax'+str(i)],Opp)
        Sub = vars()['nS'+files1[0]][Loc]
                
        n = F.nFinder(vars()['yNewMaxT'+str(i)][k],Opp,vars()['xNewMaxT'+str(i)][k],Sub)

        vars()['nForMax'+str(i)] = np.append(vars()['nForMax'+str(i)],n)
    
    for j in range(len(vars()['xNewMinT'+str(i)])):
        
        
        Loc1 = np.where(xP == F.FindNearest(xP,vars()['xNewMinT'+str(i)][j]))[0][0]
        Opp1 = vars()['yPMax'+str(i)][Loc1]   
        vars()['TMForMin'+str(i)] = np.append(vars()['TMForMin'+str(i)],Opp1)
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
            mCalc =  (2*vars()['nForMax'+str(i)][0]*vars()['dForMax'+str(i)][0]) / vars()['xNewMaxT'+str(i)][0]
            vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
        elif k == (len(vars()['xNewMaxT'+str(i)])-1):
            mCalc =  (2*vars()['nForMax'+str(i)][k]*vars()['dForMax'+str(i)][k-1]) / vars()['xNewMaxT'+str(i)][k]
            vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
        else:
            for l in range(0,2):
                mCalc =  (2*vars()['nForMax'+str(i)][k]*vars()['dForMax'+str(i)][k-1+l]) / vars()['xNewMaxT'+str(i)][k]
                vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
    
    for j in range(len(vars()['xNewMinT'+str(i)])):
        if j == 0:
            mCalc1 =  (2*vars()['nForMin'+str(i)][0]*vars()['dForMin'+str(i)][0]) / vars()['xNewMinT'+str(i)][0]
            vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))
        elif j == (len(vars()['xNewMinT'+str(i)])-1):
            mCalc1 =  (2*vars()['nForMin'+str(i)][j]*vars()['dForMin'+str(i)][j-1]) / vars()['xNewMinT'+str(i)][j]
            vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))
        else:
            for o in range(0,2):
                mCalc1 =  (2*vars()['nForMin'+str(i)][j]*vars()['dForMin'+str(i)][j-1+o]) / vars()['xNewMinT'+str(i)][j]
                vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))


    for k in range(len(vars()['xNewMaxT'+str(i)])):
        if k == 0:
            dCalc = (vars()['mForMax'+str(i)][0] * vars()['xNewMaxT'+str(i)][0]) / (2*vars()['nForMax'+str(i)][0])
            vars()['New_dForMax'+str(i)] = np.append(vars()['New_dForMax'+str(i)],dCalc)
        elif k == (len(vars()['nForMax'+str(i)])-1):
            dCalc = (vars()['mForMax'+str(i)][k+2] * vars()['xNewMaxT'+str(i)][k]) / (2*vars()['nForMax'+str(i)][k])
            vars()['New_dForMax'+str(i)] = np.append(vars()['New_dForMax'+str(i)],dCalc)
        else:
            for l in range(0,2):       
                dCalc = (vars()['mForMax'+str(i)][k+l+HasLoopMax1] * vars()['xNewMaxT'+str(i)][k]) / (2*vars()['nForMax'+str(i)][k])
                vars()['New_dForMax'+str(i)] = np.append(vars()['New_dForMax'+str(i)],dCalc)
            if l == 1: 
                HasLoopMax1 +=1
            else:
                pass 

    for j in range(len(vars()['xNewMinT'+str(i)])):
        if j == 0:
            dCalc = (vars()['mForMin'+str(i)][0] * vars()['xNewMinT'+str(i)][0]) / (2*vars()['nForMin'+str(i)][0])
            vars()['New_dForMin'+str(i)] = np.append(vars()['New_dForMin'+str(i)],dCalc)
        elif j == (len(vars()['nForMin'+str(i)])-1):
            dCalc = (vars()['mForMin'+str(i)][j+2] * vars()['xNewMinT'+str(i)][j]) / (2*vars()['nForMin'+str(i)][j])
            vars()['New_dForMin'+str(i)] = np.append(vars()['New_dForMin'+str(i)],dCalc)
        else:
            for o in range(0,2):       
                dCalc = (vars()['mForMin'+str(i)][j+o+HasLoopMin1] * vars()['xNewMinT'+str(i)][j]) / (2*vars()['nForMin'+str(i)][j])
                vars()['New_dForMin'+str(i)] = np.append(vars()['New_dForMin'+str(i)],dCalc)
            if o == 1: 
                HasLoopMin1 +=1
            else:
                pass 
 
    
    # for k in range(len(vars()['xNewMaxT'+str(i)])):
    #     if k == 0:
    #         nCalc = vars()['mForMax'+str(i)][0]*vars()['xNewMaxT'+str(i)][0] / (2*vars()['New_dForMax'+str(i)][0])
    #         vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
    #     elif k == (len(vars()['nForMax'+str(i)])-1):
    #         nCalc = vars()['mForMax'+str(i)][k+2]* vars()['xNewMaxT'+str(i)][k] / (2*vars()['New_dForMax'+str(i)][k+2])
    #         vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
    #     else:
    #         for l in range(0,2):       
    #             nCalc = vars()['mForMax'+str(i)][k+l+HasLoopMax2]* vars()['xNewMaxT'+str(i)][k] / (2*vars()['New_dForMax'+str(i)][k+l+HasLoopMax2])
    #             vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
    #         if l == 1: 
    #             HasLoopMax2 +=1
    #         else:
    #             pass 

    # for j in range(len(vars()['xNewMinT'+str(i)])):
    #     if j == 0:
    #         nCalc = vars()['mForMin'+str(i)][0]*vars()['xNewMinT'+str(i)][0] / (2*vars()['New_dForMin'+str(i)][0])
    #         vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
    #     elif j == (len(vars()['nForMin'+str(i)])-1):
    #         nCalc = vars()['mForMin'+str(i)][j+2]* vars()['xNewMinT'+str(i)][j] / (2*vars()['New_dForMin'+str(i)][j+2])
    #         vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
    #     else:
    #         for o in range(0,2):       
    #             nCalc = vars()['mForMin'+str(i)][j+o+HasLoopMin2]* vars()['xNewMinT'+str(i)][j] / (2*vars()['New_dForMin'+str(i)][j+o+HasLoopMin2])
    #             vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
    #         if o == 1: 
    #             HasLoopMin2 +=1
    #         else:
    #             pass 
    
    
    for k in range(len(vars()['xNewMaxT'+str(i)])):
        if k == 0:
            nCalc = (vars()['dForMax'+str(i)][0]*vars()['nForMax'+str(i)][0]) / (vars()['New_dForMax'+str(i)][0])
            vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
        elif k == (len(vars()['nForMax'+str(i)])-1):
            nCalc = (vars()['dForMax'+str(i)][k-1]*vars()['nForMax'+str(i)][k]) / (vars()['New_dForMax'+str(i)][k+2])
            vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
        else:
            for l in range(0,2):       
                nCalc = (vars()['dForMax'+str(i)][k-1+l]*vars()['nForMax'+str(i)][k]) / (vars()['New_dForMax'+str(i)][k+l+HasLoopMax2])
                vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
            if l == 1: 
                HasLoopMax2 +=1
            else:
                pass 

    for j in range(len(vars()['xNewMinT'+str(i)])):
        if j == 0:
            nCalc = (vars()['dForMin'+str(i)][0]*vars()['nForMin'+str(i)][0]) / (vars()['New_dForMin'+str(i)][0])
            vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
        elif j == (len(vars()['nForMin'+str(i)])-1):
            nCalc = (vars()['dForMin'+str(i)][j-1]*vars()['nForMin'+str(i)][j]) / (vars()['New_dForMin'+str(i)][j+2])
            vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
        else:
            for o in range(0,2):       
                nCalc = (vars()['dForMin'+str(i)][j-1+o]*vars()['nForMin'+str(i)][j]) / (vars()['New_dForMin'+str(i)][j+o+HasLoopMin2])
                vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
            if o == 1: 
                HasLoopMin2 +=1
            else:
                pass 
    
    
    HasLoopMax3 = 0
    HasLoopMin3 = 0
    
    
    for k in range(len(vars()['xNewMaxT'+str(i)])):
        if k == 0:
            vars()['df'+str(i)] = pd.DataFrame({'λ':[vars()['xNewMaxT'+str(i)][0]], 
                                'TM':[vars()['yNewMaxT'+str(i)][0]], 
                                'Tm':[vars()['TmForMax'+str(i)][0]], 
                                'n1':[vars()['nForMax'+str(i)][0]], 
                                'd1':[vars()['dForMax'+str(i)][0]], 
                                'm':[vars()['mForMax'+str(i)][0]], 
                                'ma':["-"] , 
                                'd2':[vars()['New_dForMax'+str(i)][0]], 
                                'd2a':["-"] , 
                                'n2':[vars()['New_nForMax'+str(i)][0]] , 
                                'n2a':["-"]})
        elif k == (len(vars()['nForMax'+str(i)])-1):
            df1 = pd.DataFrame({'λ':[vars()['xNewMaxT'+str(i)][k]], 
                                'TM':[vars()['yNewMaxT'+str(i)][k]], 
                                'Tm':[vars()['TmForMax'+str(i)][k]], 
                                'n1':[vars()['nForMax'+str(i)][k]], 
                                'd1':[vars()['dForMax'+str(i)][k-1]], 
                                'm':[vars()['mForMax'+str(i)][k+2]], 
                                'ma':["-"] , 
                                'd2':[vars()['New_dForMax'+str(i)][k+2]], 
                                'd2a':["-"] , 
                                'n2':[vars()['New_nForMax'+str(i)][k+2]] , 
                                'n2a':["-"]})
            
            vars()['df'+str(i)] = vars()['df'+str(i)].append(df1, ignore_index = True)
        else:
            df1 = pd.DataFrame({'λ':[vars()['xNewMaxT'+str(i)][k]], 
                                'TM':[vars()['yNewMaxT'+str(i)][k]], 
                                'Tm':[vars()['TmForMax'+str(i)][k]], 
                                'n1':[vars()['nForMax'+str(i)][k]], 
                                'd1':[vars()['dForMax'+str(i)][k-1]], 
                                'm':[vars()['mForMax'+str(i)][k+HasLoopMax3]], 
                                'ma':[vars()['mForMax'+str(i)][k+1+HasLoopMax3]] , 
                                'd2':[vars()['New_dForMax'+str(i)][k+HasLoopMax3]], 
                                'd2a':[vars()['New_dForMax'+str(i)][k+1+HasLoopMax3]] , 
                                'n2':[vars()['New_nForMax'+str(i)][k+HasLoopMax3]] , 
                                'n2a':[vars()['New_nForMax'+str(i)][k+1+HasLoopMax3]]})
            vars()['df'+str(i)] = vars()['df'+str(i)].append(df1, ignore_index = True) 
            HasLoopMax3 +=1


    for j in range(len(vars()['xNewMinT'+str(i)])):
        if j == 0:
            df1 = pd.DataFrame({'λ':[vars()['xNewMinT'+str(i)][0]], 
                                'TM':[vars()['TMForMin'+str(i)][0]], 
                                'Tm':[vars()['yNewMinT'+str(i)][0]], 
                                'n1':[vars()['nForMin'+str(i)][0]], 
                                'd1':[vars()['dForMin'+str(i)][0]], 
                                'm':[vars()['mForMin'+str(i)][0]],
                                'ma':["-"],
                                'd2':[vars()['New_dForMin'+str(i)][0]], 
                                'd2a':["-"] , 
                                'n2':[vars()['New_nForMin'+str(i)][0]] , 
                                'n2a':["-"]})
            vars()['df'+str(i)] = vars()['df'+str(i)].append(df1, ignore_index = True)
        elif j == (len(vars()['nForMin'+str(i)])-1):
            df1 = pd.DataFrame({'λ':[vars()['xNewMinT'+str(i)][j]], 
                                'TM':[vars()['TMForMin'+str(i)][j]], 
                                'Tm':[vars()['yNewMinT'+str(i)][j]], 
                                'n1':[vars()['nForMin'+str(i)][j]], 
                                'd1':[vars()['dForMin'+str(i)][j-1]], 
                                'm':[vars()['mForMin'+str(i)][j+2]], 
                                'ma':["-"] ,
                                'd2':[vars()['New_dForMin'+str(i)][j+2]], 
                                'd2a':["-"] , 
                                'n2':[vars()['New_nForMin'+str(i)][j+2]] , 
                                'n2a':["-"]})
            
            vars()['df'+str(i)] = vars()['df'+str(i)].append(df1, ignore_index = True)
        else:
            df1 = pd.DataFrame({'λ':[vars()['xNewMinT'+str(i)][j]], 
                                'TM':[vars()['TMForMin'+str(i)][j]], 
                                'Tm':[vars()['yNewMinT'+str(i)][j]], 
                                'n1':[vars()['nForMin'+str(i)][j]], 
                                'd1':[vars()['dForMin'+str(i)][j-1]], 
                                'm':[vars()['mForMin'+str(i)][k+HasLoopMin3]], 
                                'ma':[vars()['mForMin'+str(i)][k+1+HasLoopMin3]] , 
                                'd2':[vars()['New_dForMin'+str(i)][k+HasLoopMin3]], 
                                'd2a':[vars()['New_dForMin'+str(i)][k+1+HasLoopMin3]] , 
                                'n2':[vars()['New_nForMin'+str(i)][k+HasLoopMin3]] , 
                                'n2a':[vars()['New_nForMin'+str(i)][k+1+HasLoopMin3]]})
            vars()['df'+str(i)] = vars()['df'+str(i)].append(df1, ignore_index = True) 
            HasLoopMin3 +=1
            
    vars()['df'+str(i)] = vars()['df'+str(i)].sort_values('λ')
    # x1 = F.xFinder(n1,1.78,62.3486/100,59.6286/100)
    # x2 = F.xFinder(n2,1.775,63.6930/100,60.8116/100)

    # k1 = (-446e-9 / (4 * np.pi * d)) * np.log(x1)
    # k2 = (-534e-9 / (4 * np.pi * d)) * np.log(x2)

