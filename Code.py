import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
import csv

#dt = 1


#Importing and setting up data for processing

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))

xMin = 200
xMax = 800

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
suffix = ".txt"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)



    vars()['yFiltered'+str(i)] = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])
    
        
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)] = F.FindAntiNode(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)])
    
    vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)],vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)] = F.BackNodeFix(3,0.075,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)],vars()['MinY'+str(i)])
    
    
    #Put case statment here in python 3.10
    if i == 0:
        #All the xValues are the same, only need to do this once
        xP = np.linspace(min(vars()[files[i]+'T'][0]),max(vars()[files[i]+'T'][0]),10001)
    else:
        pass

    vars()['yPMax'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)],False)
    vars()['yPMin'+str(i)]  = F.FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)],True)


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
    # plt.figure(i)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'T'][0],vars()['yFiltered'+str(i)], label = "Filtered Data")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    
    np.savetxt(str(files[i])+".csv", vars()[files[i]+'T'].T, delimiter=',')


    plt.scatter(vars()['xNewMax'+str(i)] , vars()['yNewMax'+str(i)], color = 'black', marker = "x", label = "Maxima")
    plt.scatter(vars()['xNewMin'+str(i)] , vars()['yNewMin'+str(i)], color = 'red', marker = "x", label = "Minima")
    plt.plot(xP,vars()['yPMax'+str(i)], color = 'black', linestyle="dotted", label = "TM")
    plt.plot(xP,vars()['yPMin'+str(i)], color = 'red', linestyle="dotted", label = "Tm")
    
for i in range(len(files1)):   
    plt.plot(xP,vars()['yP'+files1[i]], color = 'orange', label = "Substrate")
    
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
    
        
    for k in range(len(vars()['xNewMax'+str(i)])):
    
        
        Loc = np.where(xP == F.FindNearest(xP,vars()['xNewMax'+str(i)][k]))[0][0]
        Opp = vars()['yPMin'+str(i)][Loc]   
        vars()['TmForMax'+str(i)] = np.append(vars()['TmForMax'+str(i)],Opp)
        Sub = vars()['nS'+files1[0]][Loc]
                
        n = F.nFinder(vars()['yNewMax'+str(i)][k],Opp,vars()['xNewMax'+str(i)][k],Sub)

        vars()['nForMax'+str(i)] = np.append(vars()['nForMax'+str(i)],n)
    

    for k in range(1,len(vars()['xNewMax'+str(i)])):
        
        d = F.dFinder(vars()['nForMax'+str(i)][k-1],vars()['nForMax'+str(i)][k],vars()['xNewMax'+str(i)][k-1],vars()['xNewMax'+str(i)][k],1)
    
        vars()['dForMax'+str(i)] = np.append(vars()['dForMax'+str(i)],abs(d))    
    
    dFin = F.dFinder(vars()['nForMax'+str(i)][0],vars()['nForMax'+str(i)][len(vars()['nForMax'+str(i)])-1],vars()['xNewMax'+str(i)][0],vars()['xNewMax'+str(i)][len(vars()['xNewMax'+str(i)])-1],len(vars()['xNewMax'+str(i)]))
    vars()['dForMax'+str(i)] = np.append(vars()['dForMax'+str(i)],abs(dFin))
    
    for k in range(len(vars()['xNewMax'+str(i)])):

        mCalc =  (2*vars()['nForMax'+str(i)][k]*vars()['dForMax'+str(i)][k]) / vars()['xNewMax'+str(i)][k]
        vars()['mForMax'+str(i)] = np.append(vars()['mForMax'+str(i)],np.around(mCalc,0))
 
        dCalc = (vars()['mForMax'+str(i)][k] * vars()['xNewMax'+str(i)][k]) / (2*vars()['nForMax'+str(i)][k])
        vars()['New_dForMax'+str(i)] = np.append(vars()['New_dForMax'+str(i)],dCalc) 
 
       
    for j in range(len(vars()['xNewMin'+str(i)])):
        
        
        Loc1 = np.where(xP == F.FindNearest(xP,vars()['xNewMin'+str(i)][j]))[0][0]
        Opp1 = vars()['yPMax'+str(i)][Loc1]   
        vars()['TMForMin'+str(i)] = np.append(vars()['TMForMin'+str(i)],Opp1)
        Sub1 = vars()['nS'+files1[0]][Loc1]
                
        n1 = F.nFinder(Opp1,vars()['yNewMin'+str(i)][j],vars()['xNewMin'+str(i)][j],Sub1)    
  
        vars()['nForMin'+str(i)] = np.append(vars()['nForMin'+str(i)],n1)
  
    
    for j in range(1,len(vars()['xNewMin'+str(i)])):
            
        d1 = F.dFinder(vars()['nForMin'+str(i)][j-1],vars()['nForMin'+str(i)][j],vars()['xNewMin'+str(i)][j-1],vars()['xNewMin'+str(i)][j],1)
    
        vars()['dForMin'+str(i)] = np.append(vars()['dForMin'+str(i)],abs(d1))
        
    dFin1 = F.dFinder(vars()['nForMin'+str(i)][0],vars()['nForMin'+str(i)][len(vars()['nForMin'+str(i)])-1],vars()['xNewMin'+str(i)][0],vars()['xNewMin'+str(i)][len(vars()['xNewMin'+str(i)])-1],len(vars()['xNewMin'+str(i)]))
    vars()['dForMin'+str(i)] = np.append(vars()['dForMin'+str(i)],abs(dFin1))    

    
    for j in range(len(vars()['xNewMin'+str(i)])):
        mCalc1 =  (2*vars()['nForMin'+str(i)][j]*vars()['dForMin'+str(i)][j]) / vars()['xNewMin'+str(i)][j]
        vars()['mForMin'+str(i)] = np.append(vars()['mForMin'+str(i)],F.round_off_rating(mCalc1))


        dCalc = (vars()['mForMin'+str(i)][j] * vars()['xNewMin'+str(i)][j]) / (2*vars()['nForMin'+str(i)][j])
        vars()['New_dForMin'+str(i)] = np.append(vars()['New_dForMin'+str(i)],dCalc)

    
    d2Combi = np.append(vars()['New_dForMax'+str(i)],vars()['New_dForMin'+str(i)])
    LambdaCombi = np.append(vars()['xNewMax'+str(i)],vars()['xNewMin'+str(i)])

    dAvg, dError, RejectedLambdas = F.ThicknessAcceptance(LambdaCombi,d2Combi, 2) 
   

    for k in range(len(vars()['xNewMax'+str(i)])):
        
        nCalc = (vars()['mForMax'+str(i)][k] * vars()['xNewMax'+str(i)][k]) / (2*dAvg)
        #nCalc = (vars()['dForMax'+str(i)][k-1+l]*vars()['nForMax'+str(i)][k]) / (vars()['New_dForMax'+str(i)][k+l+HasLoopMax2])
        vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)
    
    for j in range(len(vars()['xNewMin'+str(i)])):
        nCalc = (vars()['mForMin'+str(i)][j]*vars()['xNewMin'+str(i)][j]) / (2*dAvg)
        #  nCalc = (vars()['dForMin'+str(i)][j-1+o]*vars()['nForMin'+str(i)][j]) / (vars()['New_dForMin'+str(i)][j+o+HasLoopMin2])
        vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)
    
    
    TMcombi = np.append(vars()['yNewMax'+str(i)],vars()['TMForMin'+str(i)])
    TmCombi = np.append(vars()['TmForMax'+str(i)],vars()['yNewMin'+str(i)])
    d1Combi = np.append(vars()['dForMax'+str(i)],vars()['dForMin'+str(i)])
    n1Combi = np.append(vars()['nForMax'+str(i)],vars()['nForMin'+str(i)])
    mCombi = np.append(vars()['mForMax'+str(i)],vars()['mForMin'+str(i)])
    n2Combi = np.append(vars()['New_nForMax'+str(i)],vars()['New_nForMin'+str(i)])
    
    LambdaInds = LambdaCombi.argsort()
    
    LambdaCombi = LambdaCombi[LambdaInds[::1]]
    TMcombi = TMcombi[LambdaInds[::1]]
    TmCombi = TmCombi[LambdaInds[::1]]
    n1Combi = n1Combi[LambdaInds[::1]]
    d1Combi = d1Combi[LambdaInds[::1]]
    mCombi = mCombi[LambdaInds[::1]]
    d2Combi = d2Combi[LambdaInds[::1]]
    n2Combi = n2Combi[LambdaInds[::1]]
    
    
    header = ['λ','TM','Tm','n1','d1','m','d2','n2','discarded']
    header1 = ['dAvg','dError']
    
    with open('CSVout/'+files[i]+ '.csv','w', encoding='utf-8-sig', newline='') as f:
        
        writer = csv.writer(f)
        
        writer.writerow(header)
    
        for k in range(len(LambdaCombi)):
            if np.isin(LambdaCombi[k],RejectedLambdas) == True:
                
                data = [LambdaCombi[k],TMcombi[k],TmCombi[k], 
                       n1Combi[k],d1Combi[k],mCombi[k], 
                       d2Combi[k],n2Combi[k],"yes"]
                writer.writerow(data)

            else:
                
                data = [LambdaCombi[k],TMcombi[k],TmCombi[k], 
                       n1Combi[k],d1Combi[k],mCombi[k], 
                       d2Combi[k],n2Combi[k],"no"]
                writer.writerow(data)

    with open('CSVout/'+files[i]+ '_additional_info.csv','w', encoding='utf-8-sig', newline='') as f1:
        
        writer = csv.writer(f1)
        writer.writerow(header1)
        data = [dAvg,dError]
        writer.writerow(data)
        
        
        

            
    # vars()['df'+str(i)] = vars()['df'+str(i)].sort_values('λ')
    # x1 = F.xFinder(n1,1.78,62.3486/100,59.6286/100)
    # x2 = F.xFinder(n2,1.775,63.6930/100,60.8116/100)

    # k1 = (-446e-9 / (4 * np.pi * d)) * np.log(x1)
    # k2 = (-534e-9 / (4 * np.pi * d)) * np.log(x2)

