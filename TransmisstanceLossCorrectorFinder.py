import numpy as np
import os
import matplotlib.pyplot as plt
import Functions as F
from scipy.stats import linregress

#dt = 1


#Importing and setting up data for processing

# path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\SubstrateData\\LabData\\DP'))
# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\SubstrateData\\LabData\\Normal'))

# path2, dirs2, files2 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
suffix = ".txt"

# xMax = 800
# xMin = 200

# for i in range(len(files2)):
#     files2[i] = files2[i].rstrip(".txt")
    
#     vars()[files2[i]] = np.loadtxt(open(path2 + "\\" + files2[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
#     vars()[files2[i]][0] = vars()[files2[i]][0] * 1000
#     vars()[files2[i]][1] = vars()[files2[i]][1] / 100
    
#     vars()[files2[i]+'T'] = F.DYThorLabs(F.Trim(vars()[files2[i]],xMin,xMax))

# if len(files1) == len(files):
    
#     for i in range(len(files)):
#         files[i] = files[i][:-len(suffix)]  
        
#         vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
#         vars()[files[i]][1] = vars()[files[i]][1] / 100
        
#         files1[i] = files1[i][:-len(suffix)]  
        
#         vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
#         vars()[files1[i]][1] = vars()[files1[i]][1] / 100
        
#         Multiplier = vars()[files[i]][1] / vars()[files1[i]][1] 
        
#         np.savetxt('Multiplier.txt',Multiplier,delimiter = ',')
        
#         plt.figure(1)
#         plt.plot(vars()[files[i]][0],vars()[files[i]][1],label = 'Double Polished')
#         plt.plot(vars()[files1[i]][0],vars()[files1[i]][1],label = 'Normal')
#         plt.plot(vars()[files2[i]][0],vars()[files2[i]][1],label = 'ThorLabs')
#         plt.plot(vars()[files1[i]][0],vars()[files1[i]][1] * (max(vars()[files[i]][1])/max(vars()[files1[i]][1])),label = 'Normal Scaled')
#         plt.legend()

#         xFit = np.linspace(190,800,1001)
#         DPFit = linregress(vars()[files[i]][0][10:],vars()[files[i]][1][10:])
#         NormFit = linregress(vars()[files1[i]][0][110:],vars()[files1[i]][1][110:])
#         ThorLabsFit = linregress(vars()[files2[i]][0],vars()[files2[i]][1])
        
#         yDP = DPFit[0]*xFit + DPFit[1]
#         yNormFit = NormFit[0]*xFit + NormFit[1]
#         yThorLabs = ThorLabsFit[0]*xFit + ThorLabsFit[1]
        
#         plt.figure(2)
#         plt.plot(xFit,yDP, label = 'DP')
#         plt.plot(xFit,yNormFit, label = 'Normal')
#         plt.plot(xFit,yThorLabs, label = 'Thor Labs')
#         plt.legend()
        
#         xFit2 = np.arange(190,801,1)
#         y1 = NormFit[0]*xFit2 + NormFit[1]
#         y2 = DPFit[0]*xFit2 + DPFit[1]
        
#         Multiplier2 = y2/y1
        
#         np.savetxt('Multiplier2.txt',Multiplier2,delimiter = ',')
        
#         print('Delta is' + str(DPFit[0] - NormFit[0]))
   
# path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\DataHidden\\260822\\Standard'))
# path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\DataHidden\\260822\\V6'))


 
# if len(files1) == len(files):
#     BigArray = np.array([[0] for _ in range(len(files))])
#     for i in range(len(files)):
#         files[i] = files[i][:-len(suffix)]  
#         vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
#         vars()[files[i]][1] = vars()[files[i]][1] / 100
        
#         files1[i] = files1[i][:-len(suffix)]  
        
#         vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
#         vars()[files1[i]][1] = vars()[files1[i]][1] / 100
        
#         vars()['Multiplier'+files[i]] = vars()[files[i]][1] / vars()[files1[i]][1] 
        
#         #BigArray = np.concatenate((BigArray,vars()['Multiplier'+files[i]]),0)
        
#         np.savetxt('Multiplier.txt'+str(i),vars()['Multiplier'+files[i]],delimiter = ',')
        
#         plt.figure(1)
#         plt.plot(vars()[files[i]][0],vars()[files[i]][1],label = 'Standard')
#         plt.plot(vars()[files1[i]][0],vars()[files1[i]][1],label = 'V6')
#         plt.plot(vars()[files1[i]][0],vars()[files1[i]][1] * (max(vars()[files[i]][1])/max(vars()[files1[i]][1])),label = 'V6 Scaled')
#         plt.legend()

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\DataHidden\\060922'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('TransmisstanceLossCorrectorFinder.py')) + '\\DataHidden\\190722'))

Multi = np.loadtxt(open("MultiplierV6.txt", "rb"), delimiter=",").T    
 

for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100
        
    #BigArray = np.concatenate((BigArray,vars()['Multiplier'+files[i]]),0)
    
    # np.savetxt('Multiplier.txt'+str(i),vars()['Multiplier'+files[i]],delimiter = ',')
    
    plt.figure(1)
    plt.plot(vars()[files[i]][0],Multi*vars()[files[i]][1],label = 'Standard')
    plt.legend()