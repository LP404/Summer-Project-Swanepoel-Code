import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import csv

#dt = 1

h = 6.63e-34
c = 3e8

#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))
# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData\\LabData\\DP'))

#In nm
Slitwidth = 2

xMin = 200
xMax = 800
EnableMultiplier = 2

#Disable multiplier before using
EnableCorrectionForTesting = True

GaNenable = False

#Semillier for Sapphire
# A1 = 1.4313493
# A2 = 0.65054713
# A3 = 5.3414021
# B1 = 0.0726631
# B2 = 0.1193242
# B3 = 18.028251

A1 = 1.4313493
A2 = 0.65054713
A3 = 5.3414021
B1 = 5.2799261e-3
B2 = 1.42382647e-2
B3 = 325.017834

#Transmission unceranity
TU = (0.5/100)
#In nm
DeltaLam = 0.3

if EnableMultiplier == 0: 
    Multi = np.loadtxt(open("MultiplierV6.txt", "rb"), delimiter=",").T    
elif EnableMultiplier == 1: 
    Multi = 1.254378294491915
else:
    Multi = 1


#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
suffix = ".txt"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    vars()[files[i]][1] = Multi * vars()[files[i]][1]

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)


    #Change 0.075 back to 0.2
    vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])
    
    
    
    if i == 0:
        #All the xValues are the same, only need to do this once
        xP = np.linspace(min(vars()[files[i]+'T'][0]),max(vars()[files[i]+'T'][0]),10001)
    else:
        pass

for i in range(len(files1)):
    files1[i] = files1[i].rstrip(".txt")
    
    vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    
    #This can be commented out with an alternet approach
    vars()[files1[i]][0] = vars()[files1[i]][0]
    
    vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
    # vars()[files1[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files1[i]],xMin,xMax))
    vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],xMin,xMax)

    vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])

    vars()['nS'+files1[i]] = F.Sub_nFinder(vars()['yP'+files1[i]])
    
    if GaNenable == True:
        vars()[f"MaxValSubT_{i}"] = max(vars()['yP'+files1[i]])
    else:
        pass
    
    # xP2 = xP/1000
    
    # vars()['nS'+files1[i]] = np.sqrt((1 + ((A1*(xP2**2) / ((xP2**2) - (B1**2))) + ((A2*(xP2**2)) / ((xP2**2) - (B2**2))) + ((A3*(xP2**2)) / ((xP2**2) - (B3**2))))))

if GaNenable == True:    
    for i in range(len(files)):
        
        k = vars()[f"MaxValSubT_{0}"] / max(vars()['yFiltered'+str(i)])
        
        vars()['yFiltered'+str(i)] = k * vars()['yFiltered'+str(i)]
        
        

else:
    pass


if EnableCorrectionForTesting  == True:    
    for i in range(len(files)):
        
        MaxLocOriginal = np.where(vars()['yFiltered'+str(i)] == max(vars()['yFiltered'+str(i)]))[0][0]
        
        Val = vars()[files[0]+'T'][0][MaxLocOriginal]
        
        MaxLocSub = np.where(xP == F.FindNearestVal(xP,Val))[0][0]
        
        SubValue = vars()['yP'+files1[0]][MaxLocSub]
        
        CorrectiveMultiplier = vars()['yP'+files1[0]][MaxLocSub] / vars()['yFiltered'+str(i)][MaxLocOriginal]
        
        vars()['yFiltered'+str(i)] =  vars()['yFiltered'+str(i)] * CorrectiveMultiplier

else:
    pass

for i in range(len(files)):           
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])
    
    
    #0.075 is standard value for the second variable in the function. but adjust as neccessary
    vars()['xNewMax'+str(i)], vars()['yNewMaxUnCorr'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMinUnCorr'+str(i)] = F.BackNodeFix(3,0.075,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)])
    
    
    vars()['MaxW'+str(i)] = np.array([])
    vars()['MinW'+str(i)] = np.array([])
    

    if (max(vars()['xNewMax'+str(i)]) > max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) < min(vars()['xNewMin'+str(i)])):
        for k in range(len(vars()['xNewMin'+str(i)])):
            if k == 0:
                W = 2 * (min(vars()['xNewMin'+str(i)]) - min(vars()['xNewMax'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)
    
        W = 2 * (max(vars()['xNewMax'+str(i)]) - max(vars()['xNewMin'+str(i)]))   
        vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)
        
        for j in range(1,len(vars()['xNewMax'+str(i)])):
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w) 
                
    elif (max(vars()['xNewMax'+str(i)]) < max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) < min(vars()['xNewMin'+str(i)])):
        for k in range(len(vars()['xNewMin'+str(i)])):
            if k == 0:
                W = 2 * (min(vars()['xNewMin'+str(i)]) - min(vars()['xNewMax'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)      
                
        for j in range(1,len(vars()['xNewMax'+str(i)])):
            if j == (len(vars()['xNewMax'+str(i)]) - 1):        
                w = 2 * (max(vars()['xNewMin'+str(i)]) - max(vars()['xNewMax'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
                
    elif (max(vars()['xNewMax'+str(i)]) > max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) > min(vars()['xNewMin'+str(i)])):
        for k in range(1,len(vars()['xNewMin'+str(i)])):
            if k == (len(vars()['xNewMin'+str(i)]) - 1):
                W = 2 * (max(vars()['xNewMax'+str(i)]) - max(vars()['xNewMin'+str(i)]))
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)              
            else:
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
        
        for j in range(1,len(vars()['xNewMax'+str(i)])):
            if j == (len(vars()['xNewMax'+str(i)]) - 1):
                w = 2 * (min(vars()['xNewMax'+str(i)]) - min(vars()['xNewMin'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)   
            else:
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
    
    elif (max(vars()['xNewMax'+str(i)]) < max(vars()['xNewMin'+str(i)])) and (min(vars()['xNewMax'+str(i)]) > min(vars()['xNewMin'+str(i)])):
        for k in range(1,len(vars()['xNewMin'+str(i)])):
                W = vars()['xNewMin'+str(i)][k] - vars()['xNewMin'+str(i)][k-1]
                vars()['MaxW'+str(i)] = np.append(vars()['MaxW'+str(i)],W)           
        
        for j in range(len(vars()['xNewMax'+str(i)])):
            if j == 0:
                w = 2 * (min(vars()['xNewMax'+str(i)]) - min(vars()['xNewMin'+str(i)]))
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)   
            else:
                w = vars()['xNewMax'+str(i)][j] - vars()['xNewMax'+str(i)][j-1]
                vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w) 
            
        w = 2 * (max(vars()['xNewMin'+str(i)]) - max(vars()['xNewMax'+str(i)]))   
        vars()['MinW'+str(i)] = np.append(vars()['MinW'+str(i)],w)
        
    
    vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)] + ((Slitwidth*vars()['yNewMaxUnCorr'+str(i)])/vars()['MaxW'+str(i)])**2
    vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)] + ((Slitwidth*vars()['yNewMinUnCorr'+str(i)])/vars()['MinW'+str(i)])**2
 

    # vars()['yNewMax'+str(i)] = vars()['yNewMaxUnCorr'+str(i)]
    # vars()['yNewMin'+str(i)] = vars()['yNewMinUnCorr'+str(i)]


    # np.savetxt('xP.txt', xP, delimiter=",")
    # np.savetxt('x.txt', vars()[files[i]+'T'][0], delimiter=",")
    # np.savetxt('y.txt',vars()['yFiltered'+str(i)], delimiter=",")
    # np.savetxt('AntiX.txt',vars()['xNewMax'+str(i)], delimiter=",")
    # np.savetxt('AntiY.txt',vars()['yNewMaxUnCorr'+str(i)], delimiter=",")
    
    vars()['yPMaxUnCorr'+str(i)]  = F.Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMaxUnCorr'+str(i)])
    vars()['yPMinUnCorr'+str(i)]  = F.Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMinUnCorr'+str(i)])

    # vars()['yPMaxUnCorr'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMaxUnCorr'+str(i)],False)
    # vars()['yPMinUnCorr'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMaxUnCorr'+str(i)],True)

    # vars()['yPMax'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)],False)
    # vars()['yPMin'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)],True)


    vars()['yPMax'+str(i)] = F.Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)])
    vars()['yPMin'+str(i)] = F.Reinterp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)])

for i in range(len(files)):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.figure(i,figsize = (16,10),dpi = 600)
    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)], label = "Filtered Data")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    
    # np.savetxt(str(files[i])+".csv", vars()[files[i]+'T'].T, delimiter=',')


    plt.scatter(vars()['xNewMax'+str(i)] , vars()['yNewMaxUnCorr'+str(i)], color = 'black', marker = "x", label = "Maxima")
    plt.scatter(vars()['xNewMin'+str(i)] , vars()['yNewMinUnCorr'+str(i)], color = 'red', marker = "x", label = "Minima")
    plt.plot(xP,vars()['yPMaxUnCorr'+str(i)], color = 'black', linestyle="dotted", label = "TM")
    plt.plot(xP,vars()['yPMinUnCorr'+str(i)], color = 'red', linestyle="dotted", label = "Tm")
    plt.plot(xP,vars()['yP'+files1[0]], color = 'orange', label = "Substrate")
    
    plt.legend()
    
for i in range(len(files)):
    

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

    vars()['TmForMaxUnCorr'+str(i)] = np.array([])
    vars()['TMForMinUnCorr'+str(i)] = np.array([])
    
    vars()['n1UncertMax'+str(i)] = np.array([])
    vars()['n1UncertMin'+str(i)] = np.array([])
    
    vars()['d1UncertMax'+str(i)] = np.array([])
    vars()['d1UncertMin'+str(i)] = np.array([]) #f'file_{i}'
       
    
    #Finds values that do not account for slit width correction
    #These values are not used for any calcultions beyond finding corrected values
    
    #dUncert(Lam1,Lam2,n1,n2,d,LamUncert,n1Uncert,n2Uncert)
    #nUncertFinder(TUncert,TM,Tm,n)
    
    for k in range(len(vars()['xNewMax'+str(i)])):
    
        
        Loc = np.where(xP == F.FindNearestVal(xP,vars()['xNewMax'+str(i)][k]))[0][0]
        Opp = vars()['yPMinUnCorr'+str(i)][Loc]   
        vars()['TmForMaxUnCorr'+str(i)] = np.append(vars()['TmForMaxUnCorr'+str(i)],Opp)

    for j in range(len(vars()['xNewMin'+str(i)])):
        
        
        Loc1 = np.where(xP == F.FindNearestVal(xP,vars()['xNewMin'+str(i)][j]))[0][0]
        Opp1 = vars()['yPMaxUnCorr'+str(i)][Loc1]   
        vars()['TMForMinUnCorr'+str(i)] = np.append(vars()['TMForMinUnCorr'+str(i)],Opp1)
        
        
    TMcombiUnCorr = np.append(vars()['yNewMaxUnCorr'+str(i)],vars()['TMForMinUnCorr'+str(i)])
    TmCombiUnCorr = np.append(vars()['TmForMaxUnCorr'+str(i)],vars()['yNewMinUnCorr'+str(i)])
    
        
    for k in range(len(vars()['xNewMax'+str(i)])):
    
        
        Loc2 = np.where(xP == F.FindNearestVal(xP,vars()['xNewMax'+str(i)][k]))[0][0]
        Opp2 = vars()['yPMin'+str(i)][Loc2]   
        vars()['TmForMax'+str(i)] = np.append(vars()['TmForMax'+str(i)],Opp2)
        Sub = vars()['nS'+files1[0]][Loc2]
                
        n = F.nFinder(vars()['yNewMax'+str(i)][k],Opp2,vars()['xNewMax'+str(i)][k],Sub)
 
        vars()['nForMax'+str(i)] = np.append(vars()['nForMax'+str(i)],n)

        nError = F.nUncertFinder(TU,vars()['xNewMax'+str(i)][k],Opp2,n)
        
        vars()['n1UncertMax'+str(i)] = np.append(vars()['n1UncertMax'+str(i)],nError)

    for k in range(1,len(vars()['xNewMax'+str(i)])):
        
        d = F.dFinder(vars()['nForMax'+str(i)][k-1],vars()['nForMax'+str(i)][k],vars()['xNewMax'+str(i)][k-1],vars()['xNewMax'+str(i)][k],1)
    
        vars()['dForMax'+str(i)] = np.append(vars()['dForMax'+str(i)],abs(d))    
    
        #dErr = F.dUncert(vars()['xNewMax'+str(i)][k-1],vars()['xNewMax'+str(i)][k],vars()['nForMax'+str(i)][k-1],vars()['nForMax'+str(i)][k],d,DeltaLam,vars()['n1UncertMax'+str(i)][k-1],vars()['n1UncertMax'+str(i)][k])
   
        dErr = F.dUncerIgnN(vars()['xNewMax'+str(i)][k-1],vars()['xNewMax'+str(i)][k],vars()['nForMax'+str(i)][k-1],vars()['nForMax'+str(i)][k],d,DeltaLam)
    
        vars()['d1UncertMax'+str(i)] = np.append(vars()['d1UncertMax'+str(i)],dErr)
   
    
    dFin = F.dFinder(vars()['nForMax'+str(i)][0],vars()['nForMax'+str(i)][len(vars()['nForMax'+str(i)])-1],vars()['xNewMax'+str(i)][0],vars()['xNewMax'+str(i)][len(vars()['xNewMax'+str(i)])-1],len(vars()['xNewMax'+str(i)]-1))
    vars()['dForMax'+str(i)] = np.append(vars()['dForMax'+str(i)],abs(dFin))
        
    #dErrFin = F.dUncert(vars()['xNewMax'+str(i)][0],vars()['xNewMax'+str(i)][len(vars()['xNewMax'+str(i)])-1],vars()['nForMax'+str(i)][0],vars()['nForMax'+str(i)][len(vars()['nForMax'+str(i)])-1],dFin,DeltaLam,vars()['n1UncertMax'+str(i)][0],vars()['n1UncertMax'+str(i)][len(vars()['n1UncertMax'+str(i)])-1])

    dErrFin = F.dUncerIgnN(vars()['xNewMax'+str(i)][0],vars()['xNewMax'+str(i)][len(vars()['xNewMax'+str(i)])-1],vars()['nForMax'+str(i)][0],vars()['nForMax'+str(i)][len(vars()['nForMax'+str(i)])-1],dFin,DeltaLam)

    vars()['d1UncertMax'+str(i)] = np.append(vars()['d1UncertMax'+str(i)],dErrFin)
    
    
    for j in range(len(vars()['xNewMin'+str(i)])):
        
        
        Loc3 = np.where(xP == F.FindNearestVal(xP,vars()['xNewMin'+str(i)][j]))[0][0]
        Opp3 = vars()['yPMax'+str(i)][Loc3]   
        vars()['TMForMin'+str(i)] = np.append(vars()['TMForMin'+str(i)],Opp3)
        Sub1 = vars()['nS'+files1[0]][Loc3]
                
        n1 = F.nFinder(Opp3,vars()['yNewMin'+str(i)][j],vars()['xNewMin'+str(i)][j],Sub1)    
  
        vars()['nForMin'+str(i)] = np.append(vars()['nForMin'+str(i)],n1)
  
        n1Error = F.nUncertFinder(TU,Opp3,vars()['xNewMin'+str(i)][j],n1)

        vars()['n1UncertMin'+str(i)] = np.append(vars()['n1UncertMin'+str(i)],n1Error)
    
    n1Combi = np.append(vars()['nForMax'+str(i)],vars()['nForMin'+str(i)])
    n1ErrCombi = np.append(vars()['n1UncertMax'+str(i)],vars()['n1UncertMin'+str(i)])
    
    for j in range(1,len(vars()['xNewMin'+str(i)])):
            
        d1 = F.dFinder(vars()['nForMin'+str(i)][j-1],vars()['nForMin'+str(i)][j],vars()['xNewMin'+str(i)][j-1],vars()['xNewMin'+str(i)][j],1)
    
        vars()['dForMin'+str(i)] = np.append(vars()['dForMin'+str(i)],abs(d1))
        
        #d1Err = F.dUncert(vars()['xNewMin'+str(i)][j-1],vars()['xNewMin'+str(i)][j],vars()['nForMin'+str(i)][j-1],vars()['nForMin'+str(i)][j],d1,DeltaLam,vars()['n1UncertMin'+str(i)][j-1],vars()['n1UncertMin'+str(i)][j])
        
        d1Err = F.dUncerIgnN(vars()['xNewMin'+str(i)][j-1],vars()['xNewMin'+str(i)][j],vars()['nForMin'+str(i)][j-1],vars()['nForMin'+str(i)][j],d1,DeltaLam)
        
        vars()['d1UncertMin'+str(i)] = np.append(vars()['d1UncertMin'+str(i)],d1Err)
        
    dFin1 = F.dFinder(vars()['nForMin'+str(i)][0],vars()['nForMin'+str(i)][len(vars()['nForMin'+str(i)])-1],vars()['xNewMin'+str(i)][0],vars()['xNewMin'+str(i)][len(vars()['xNewMin'+str(i)])-1],len(vars()['xNewMin'+str(i)]-1))
    vars()['dForMin'+str(i)] = np.append(vars()['dForMin'+str(i)],abs(dFin1))    

    #dErrFin1 = F.dUncert(vars()['xNewMin'+str(i)][0],vars()['xNewMin'+str(i)][len(vars()['xNewMin'+str(i)])-1],vars()['nForMin'+str(i)][0],vars()['nForMin'+str(i)][len(vars()['nForMin'+str(i)])-1],dFin1,DeltaLam,vars()['n1UncertMin'+str(i)][0],vars()['n1UncertMin'+str(i)][len(vars()['n1UncertMin'+str(i)])-1])
    
    dErrFin1 = F.dUncerIgnN(vars()['xNewMin'+str(i)][0],vars()['xNewMin'+str(i)][len(vars()['xNewMin'+str(i)])-1],vars()['nForMin'+str(i)][0],vars()['nForMin'+str(i)][len(vars()['nForMin'+str(i)])-1],dFin1,DeltaLam)
    
    vars()['d1UncertMin'+str(i)] = np.append(vars()['d1UncertMin'+str(i)],dErrFin1)

    
    d1Combi = np.append(vars()['dForMax'+str(i)],vars()['dForMin'+str(i)])
    LambdaCombi = np.append(vars()['xNewMax'+str(i)],vars()['xNewMin'+str(i)])
    d1ErrCombi = np.append(vars()['d1UncertMax'+str(i)],vars()['d1UncertMin'+str(i)])
    
    d1Avg, d1StandardError, RejectedLambdas1, d1ExpError = F.ThicknessAcceptance(LambdaCombi,d1Combi,1,d1ErrCombi,2) 
    d1Error = np.sqrt((d1StandardError**2 + d1ExpError**2))
    
    for k in range(len(vars()['xNewMax'+str(i)])):

        mCalc =  (2*vars()['nForMax'+str(i)][k]*d1Avg) / vars()['xNewMax'+str(i)][k]
        vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
    
        
        dCalc = (vars()['mForMax'+str(i)][k] * vars()['xNewMax'+str(i)][k]) / (2*vars()['nForMax'+str(i)][k])
        vars()['New_dForMax'+str(i)] = np.append(vars()['New_dForMax'+str(i)],dCalc) 
    

    for j in range(len(vars()['xNewMin'+str(i)])):
        mCalc1 =  (2*vars()['nForMin'+str(i)][j]*d1Avg) / vars()['xNewMin'+str(i)][j]
        vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))

        dCalc = (vars()['mForMin'+str(i)][j] * vars()['xNewMin'+str(i)][j]) / (2*vars()['nForMin'+str(i)][j])
        
        vars()['New_dForMin'+str(i)] = np.append(vars()['New_dForMin'+str(i)],dCalc)

    d2Combi = np.append(vars()['New_dForMax'+str(i)],vars()['New_dForMin'+str(i)])
    
    Per1 = F.AbsToPer(n1ErrCombi,n1Combi)
    Per2 = F.AbsToPer(d1ErrCombi,d1Combi)
    
    d2ErrCombi = np.array([])
    
    for k in range(len(LambdaCombi)):
        if Per1[k] < Per2[k]:
            d2ErrCombi = np.append(d2ErrCombi,F.PerToAbs(Per2[k],d2Combi[k]))
        else:
            d2ErrCombi = np.append(d2ErrCombi,F.PerToAbs(Per1[k],d2Combi[k]))
    
    
    d2Avg, d2StandardError, RejectedLambdasD2, d2ExpError = F.ThicknessAcceptance(LambdaCombi,d2Combi, 0,d2ErrCombi,3) 

    d2Error = np.sqrt((d2StandardError**2 + d2ExpError**2))
    
    RejectedLambdasD2 = np.append(RejectedLambdasD2,RejectedLambdas1)

    for k in range(len(vars()['xNewMax'+str(i)])):
        
        nCalc = (vars()['mForMax'+str(i)][k] * vars()['xNewMax'+str(i)][k]) / (2*d2Avg)
        #nCalc = (vars()['dForMax'+str(i)][k-1+l]*vars()['nForMax'+str(i)][k]) / (vars()['New_dForMax'+str(i)][k+l+HasLoopMax2])
        vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
    
    
    for j in range(len(vars()['xNewMin'+str(i)])):
        nCalc = (vars()['mForMin'+str(i)][j]*vars()['xNewMin'+str(i)][j]) / (2*d2Avg)
        #  nCalc = (vars()['dForMin'+str(i)][j-1+o]*vars()['nForMin'+str(i)][j]) / (vars()['New_dForMin'+str(i)][j+o+HasLoopMin2])
        vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
        
   
    TMcombi = np.append(vars()['yNewMax'+str(i)],vars()['TMForMin'+str(i)])
    TmCombi = np.append(vars()['TmForMax'+str(i)],vars()['yNewMin'+str(i)])
    mCombi = np.append(vars()['mForMax'+str(i)],vars()['mForMin'+str(i)])
    n2Combi = np.append(vars()['New_nForMax'+str(i)],vars()['New_nForMin'+str(i)])
            
    LambdaInds = LambdaCombi.argsort()
    
    LambdaCombi = LambdaCombi[LambdaInds[::1]]
    TMcombiUnCorr = TMcombiUnCorr[LambdaInds[::1]]
    TmCombiUnCorr = TmCombiUnCorr[LambdaInds[::1]]
    TMcombi = TMcombi[LambdaInds[::1]]
    TmCombi = TmCombi[LambdaInds[::1]]
    n1Combi = n1Combi[LambdaInds[::1]]
    d1Combi = d1Combi[LambdaInds[::1]]
    mCombi = mCombi[LambdaInds[::1]]
    d2Combi = d2Combi[LambdaInds[::1]]
    n2Combi = n2Combi[LambdaInds[::1]]
    
    d1ErrCombi = d1ErrCombi[LambdaInds[::1]]
    n1ErrCombi = n1ErrCombi[LambdaInds[::1]]
    
    d2ErrCombi = d2ErrCombi[LambdaInds[::1]]
    
    mErrCombi = np.array([])
    n2ErrCombi = np.array([])
    
    Per1 = F.AbsToPer(n1ErrCombi,n1Combi)
    Per2 = F.AbsToPer(d1ErrCombi,d1Combi)
    
    for k in range(len(LambdaCombi)):
        if Per1[k] < Per2[k]:
            mErrCombi = np.append(mErrCombi,F.PerToAbs(Per2[k],mCombi[k]))
            n2ErrCombi = np.append(n2ErrCombi,F.PerToAbs(Per2[k],n2Combi[k]))
        else:
            mErrCombi = np.append(mErrCombi,F.PerToAbs(Per1[k],mCombi[k]))
            n2ErrCombi = np.append(n2ErrCombi,F.PerToAbs(Per1[k],n2Combi[k]))
    
    
    Index = np.array([])
    
    for k in range(len(LambdaCombi)):
        IndexVal = np.where(LambdaCombi[k] == np.around(xP,1))[0][0]
        Index = np.append(Index,IndexVal)
        
        
    Fval = 4 * n2Combi**2 * vars()['nS'+files1[0]][Index.astype(int)] * ((TMcombi + TmCombi)/(TMcombi*TmCombi))
    x = (Fval - np.sqrt(Fval**2 - ((n2Combi**2 - 1)**3  * (n2Combi**2 - vars()['nS'+files1[0]][Index.astype(int)]**4 ))))/ ((n2Combi-1)**3 * (n2Combi-vars()['nS'+files1[0]][Index.astype(int)]**2))
    
    kVal = (-LambdaCombi / (4 * np.pi * d2Combi)) * np.log(x)
    alphaVal = ((4 * np.pi * kVal) / LambdaCombi) * 1e7
    
    alphaVal2 = -1 * (np.log(x) / d2Combi) * 1e7
    
    
    header = ['λ (nm)','TM','Tm','TM Corrected','Tm Corrected','n1','Δn1','d1 (nm)','Δd1 (nm)','m','Δm','d2 (nm)','Δd2 (nm)','n2','Δn2','κ','α1 (cm^-1)','α2 (cm^-1)','discarded']
    header1 = ['d1Avg (nm)','Δd1(nm)','d2Avg (nm)','Δd2 (nm)']
    header2 = ['λ (nm)','Traw','TM','Tm','T true','hv','n true','α1 (cm^-1)','α2 (cm^-1)','α3 (cm^-1)','k']
    
    with open('CSVout/'+files[i]+ '.csv','w', encoding='utf-8-sig', newline='') as f:
        
        writer = csv.writer(f)
        
        writer.writerow(header)
    
        for k in range(len(LambdaCombi)):
            if np.isin(LambdaCombi[k],RejectedLambdasD2) == True:
                
                data = [LambdaCombi[k],TMcombiUnCorr[k],TmCombiUnCorr[k],
                        TMcombi[k],TmCombi[k], 
                        n1Combi[k],n1ErrCombi[k],d1Combi[k],d1ErrCombi[k],mCombi[k],mErrCombi[k], 
                        d2Combi[k],d2ErrCombi[k],n2Combi[k],n2ErrCombi[k],kVal[k],alphaVal[k],alphaVal2[k],"yes"]
                
                writer.writerow(data)

            else:
                
                data = [LambdaCombi[k],TMcombiUnCorr[k],TmCombiUnCorr[k],
                        TMcombi[k],TmCombi[k], 
                        n1Combi[k],n1ErrCombi[k],d1Combi[k],d1ErrCombi[k],mCombi[k],mErrCombi[k],
                        d2Combi[k],d2ErrCombi[k],n2Combi[k],n2ErrCombi[k],kVal[k],alphaVal[k],alphaVal2[k],"no"]
                
                writer.writerow(data)

    with open('CSVout/Thickness/'+files[i]+ '_Thickness.csv','w', encoding='utf-8-sig', newline='') as f1:
        
        writer = csv.writer(f1)
        writer.writerow(header1)
        data = [d1Avg,d1Error,d2Avg,d2Error]
        writer.writerow(data)
        
    Trn = np.sqrt((vars()['yPMax'+str(i)]  * vars()['yPMin'+str(i)]))
    

    H = ((4 * vars()['nS'+files1[0]]**2) / ((vars()['nS'+files1[0]]**2 + 1) * (Trn**2))) - (((vars()['nS'+files1[0]]**2) + 1) / 2)

    
    nTrn = np.sqrt((H + np.sqrt(((H**2) - (vars()['nS'+files1[0]]**2))))) 
    
    # #hv is identical for all measuremnts 
    hv = ((h * c) / (xP * 1e-9)) * 6.242e18
    alpha3 = ((np.log((1/Trn))) / d2Avg) * 1e7
    
    k = (alpha3*xP / 4*np.pi) / 1e7 
      
    yRawInterp = np.interp(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])
    
    alphaValInterp = np.interp(xP,LambdaCombi,alphaVal)
    alphaVal2Interp = np.interp(xP,LambdaCombi,alphaVal2)
    
    with open('CSVout/Absorption/'+files[i]+ '_Absorption.csv','w', encoding='utf-8-sig', newline='') as f2:
         
        writer = csv.writer(f2)
        writer.writerow(header2)
        
        for l in range(len(Trn)):
            data = [xP[l],yRawInterp[l],vars()['yPMax'+str(i)][l],vars()['yPMin'+str(i)][l],Trn[l],hv[l],nTrn[l],alphaValInterp[l],alphaVal2Interp[l],alpha3[l],k[l]]
            writer.writerow(data)

            
    # vars()['df'+str(i)] = vars()['df'+str(i)].sort_values('λ')
    # x1 = F.xFinder(n1,1.78,62.3486/100,59.6286/100)
    # x2 = F.xFinder(n2,1.775,63.6930/100,60.8116/100)

    # k1 = (-446e-9 / (4 * np.pi * d)) * np.log(x1)
    # k2 = (-534e-9 / (4 * np.pi * d)) * np.log(x2)