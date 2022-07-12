import numpy as np
import Functions as F
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter
import os
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

#Factor in substrate
#Add in reflection?

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))
xMin = 200
xMax = 800

def BlockFinder(xArray):
    gap = np.diff(xArray)
    BlockLength = 0
    Blocks = np.array([])
    for i in range(len(gap)):
        
        if gap[i] < np.mean(gap):
            BlockLength += 1
            if i == (len(gap) - 1):
                Blocks = np.append(Blocks, (BlockLength + 1))
            else:
                pass
        else:
            Blocks = np.append(Blocks, (BlockLength + 1))
            BlockLength = 0
    
    return Blocks

#FixThisLater
def ArrayInterlace(Array1, Array2):
    Combi = np.zeros(Array1.size + Array2.size)
    
    if len(Array1) < len(Array2):
        CritLength = len(Array1)
        CritArray = Array2    
        for i, (item1, item2) in enumerate(zip(Array1, Array2)):
            Combi[i*2] = item1
            Combi[i*2+1] = item2
    
        DeleteIndis = np.arange(2*CritLength-1,len(Combi),1)
        Combi = np.delete(Combi,DeleteIndis)
        Combi = np.append(Combi,CritArray[CritLength-1:])
        
    elif len(Array1) > len(Array2):
        CritLength = len(Array2)
        CritArray = Array1
    
        for i, (item1, item2) in enumerate(zip(Array1, Array2)):
            Combi[i*2] = item1
            Combi[i*2+1] = item2
    
        DeleteIndis = np.arange(2*CritLength,len(Combi),1)
        Combi = np.delete(Combi,DeleteIndis)
        Combi = np.append(Combi,CritArray[CritLength:])
    
    else:
    
        for i, (item1, item2) in enumerate(zip(Array1, Array2)):
            Combi[i*2] = item1
            Combi[i*2+1] = item2

    return Combi

#!!!Once case statments are released, change this
#!!! Add a MaxMin Detector After 500nm
def AntiNodeHighlander(xArray,yArray,xFull,yFull):
    
    blocks = BlockFinder(xArray)
    
    xRevAntiNode = np.array([])
    yRevAntiNode = np.array([])
    
    
    for i in range(len(blocks)+1):
        if i == 0:
            pass
        else:
            start = int(sum(blocks[0:i-1]))
            end = int(sum(blocks[0:i]))
            xTrial = xArray[start:end]
            yTrial = yArray[start:end] 
            
            if max(yTrial) < 0.2:
                pass
            else:
                #!!!Find a more resource efficient way to do this
                #!!Sort the potential crankMax loop
                if len(xTrial) == 1:
                    xRevAntiNode = np.append(xRevAntiNode,xTrial)
                    yRevAntiNode = np.append(yRevAntiNode,yTrial)
                
                elif len(xTrial) > 1:
                    MinVal = 100
                    xAccept = xTrial[0]
                    yAccept = yTrial[0]
                    for k  in range(len(yTrial)):
                        loc = np.where(yFull == yTrial[k])[0][0]
                        m2ndDer = (yFull[loc+1] - (2 * yFull[loc]) + yFull[loc-1])
                        if abs(m2ndDer) < MinVal:
                            MinVal = abs(m2ndDer)
                            yAccept = yTrial[k]
                            xAccept = xTrial[k]
                        else:
                            pass
                    xRevAntiNode = np.append(xRevAntiNode,xAccept)
                    yRevAntiNode = np.append(yRevAntiNode,yAccept)
                     
                
                else:
                    print('What?')

    return xRevAntiNode, yRevAntiNode


def AntiNodeSanityChecker(xMax,yMax,xMin,yMin):
    
    MaxLoc = np.array([])
    MinLoc = np.array([])
    
    for i in range(len(yMax)):
        if yMax[i] < 0.2:
            MaxLoc = np.append(MaxLoc,i)
            
        elif xMax[i] < 100:
            MaxLoc = np.append(MaxLoc,i)
          
    for i in range(len(yMin)):
        if yMin[i] < 0.2:
            MinLoc = np.append(MinLoc,i)
   
        elif xMin[i] < 100:
            MinLoc = np.append(MinLoc,i)
 
        
    MaxLoc = np.sort(MaxLoc).astype('int32')
    MinLoc = np.sort(MinLoc).astype('int32')
    
    yMax = np.delete(yMax,MaxLoc)
    xMax = np.delete(xMax,MaxLoc)   
    yMin = np.delete(yMin,MinLoc)
    xMin = np.delete(xMin,MinLoc) 
    

        
    # if min(xMax) <= min(xMin):
    #     Combi = ArrayInterlace(xMax, (-1 * xMin))     
    # else:
    #     Combi = ArrayInterlace((-1* xMin), * xMax)
        
    # for i in range(1,len(Combi)):
    #     if (Combi[i-1]/Combi[i-1]) != (Combi[i]/Combi[i]):
    #         pass
    #     elif (Combi[i-1]/abs(Combi[i-1])) == (Combi[i]/abs(Combi[i])) == 1:  
    #         Loc = np.where(Combi[i-1] == xMax)[0][0]
    #         Loc1 = np.where(Combi[i] == xMax)[0][0]
            
    #         if yMax[Loc] <= yMax[Loc1]:
    #             yMax = np.delete(yMax,Loc)
    #             xMax = np.delete(yMax,Loc)
    #         else:
    #             yMax = np.delete(yMax,Loc1)
    #             xMax = np.delete(yMax,Loc1)   
                  
    #     elif (Combi[i-1]/abs(Combi[i-1])) == (Combi[i]/abs(Combi[i])) == -1:
    #         Loc = np.where(abs(Combi[i-1]) == xMin)[0][0]
    #         Loc1 = np.where(abs(Combi[i]) == xMin)[0][0]
    #         if yMin[Loc] >= yMin[Loc1]:                    
    #             yMin = np.delete(yMin,Loc)
    #             xMin = np.delete(yMin,Loc)
    #         else:
    #             yMin = np.delete(yMin,Loc1)
    #             xMin = np.delete(yMin ,Loc1) 
    #     else:
    #         pass
            
        # LastLoc = np.where(abs(Combi) < 500)[0][-1]
        
        # if (Combi[LastLoc]/abs(Combi[LastLoc])) == 1:
        #     IsLastMax = True
        # else:
        #     IsLastMax = False
        
        
        # for i in range(LastLoc+1,len(Combi)):
        #         if abs(Combi[i] - Combi[i-1]) < 10:
        #             if  IsLastMax == True
                    
                
                    
    return xMax, yMax, xMin, yMin    

#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i].rstrip(".txt")
    
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    vars()[files[i]][1] = vars()[files[i]][1] / 100

    #T stands for Truncated
    vars()[files[i]+'T'] = F.Trim(vars()[files[i]],xMin,xMax)

    
    y = F.NoiseFilter(3,0.2,vars()[files[i]+'T'][1])
    
    
    plt.figure(figsize=(10,5))
    plt.plot(vars()[files[i]+'T'][0], y, 'k', linewidth=1.75, label = 'cleaner signal')
    plt.legend(loc='best')
    plt.grid(True)
    
    Xmaxima = vars()[files[i]+'T'][0][argrelextrema(y, np.greater)[0]]
    Xminima = vars()[files[i]+'T'][0][argrelextrema(y, np.less)[0]]
    
    Ymaxima = y[argrelextrema(y, np.greater)[0]]
    Yminima = y[argrelextrema(y, np.less)[0]]


    plt.scatter(Xmaxima,Ymaxima, color = 'black', marker = "x")
    plt.scatter(Xminima,Yminima, color = 'red', marker = "x")
    
#Put code in here you want to hide in the spyder ide with the close/open feature

def HideCodeInSpyder():

    #Mins
    
        # for k in range(len(vars()['nForMin'+str(i)])):
        #     if k == 0:
        #         dCalc =  (2*vars()['nForMin'+str(i)][0]) / (vars()['mForMin'+str(i)][0] * vars()['NewXForMin'+str(i)][0])
        #         vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
        #     elif k == (len(vars()['nForMin'+str(i)])-1):
        #         dCalc =  (2*vars()['nForMin'+str(i)][k]) / (vars()['mForMin'+str(i)][k] * vars()['NewXForMin'+str(i)][k])
        #         vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
        #     else:
        #         for l in range(0,2):
        #             lamCalc =  2*vars()['nForMin'+str(i)][k]* vars()['dForMin'+str(i)][k-1+l] / vars()['mForMin'+str(i)][[k+l+HasLoopMin]]
        #             vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
                
        #         if l == 1: 
        #             HasLoopMin +=1
        #         else:
        #             pass  
    
    
    # m = np.array([])
    
    # lambda1 = np.array([100,200,300,400])
    # n = np.array([2.0,1.99,1.98,1.977])
    # d = np.array([300,400,500])
    
    # HasLoop = 0
    # for i in range(len(lambda1)):
    #     if i == 0:
    #         mCalc =  2*n[0]*d[0] / lambda1[0]
    #         m = np.append(m,np.around(mCalc,0))
    #     elif i == (len(lambda1)-1):
    #         mCalc =  2*n[i]*d[i-1] / lambda1[i]
    #         m = np.append(m,np.around(mCalc,0))
    #     else:
    #         for k  in range(0,2):
    #             mCalc =  2*n[i]*d[i-1+k] / lambda1[i]
    #             m = np.append(m,np.around(mCalc,0))
    #             if k == 1: 
    #                 HasLoop +=1
    #             else:
    #                 pass
    #         pass
            
    
    
    # new_d = np.array([])
    # HasLoop1 = 0
    # for i in range(len(n)):
    #     if i == 0:
    #         dCalc =  (m[0] * lambda1[0]) / (2*n[0])
    #         new_d = np.append(new_d,dCalc)
    #     elif i == (len(n)-1):
    #         dCalc =  (m[i+2] * lambda1[i]) / (2*n[i])
    #         new_d = np.append(new_d,dCalc)
    #     else:
    #         for k  in range(0,2):
    #             dCalc =  (m[i+k+HasLoop1] * lambda1[i]) / (2*n[i])
    #             print(i)
    #             print(i-1+k)
    #             print(i+k+HasLoop1)
    #             print('')
    #             new_d = np.append(new_d,dCalc)
    #             if k == 1: 
    #                 HasLoop1 +=1
    #             else:
    #                 pass
        
    # new_n = np.array([])
    # HasLoop2 = 0
    # for i in range(len(n)):
    #     if i == 0:
    #         nCalc =  (2*n[0]*d[0]) / (2*new_d[0])
    #         new_n = np.append(new_n,nCalc)
    #     elif i == (len(n)-1):
    #         nCalc =  (2*n[i]*d[i-1]) / (2*new_d[i+2])
    #         new_n = np.append(new_n,nCalc)
    #     else:
    #         for k  in range(0,2):
    #             nCalc =  (2*n[i]*d[i-1+k]) / (2*new_d[i+k+HasLoop2])
    #             print(i)
    #             print(i-1+k)
    #             print(i+k+HasLoop1)
    #             print('')
    #             new_n = np.append(new_n,nCalc)
    #             if k == 1: 
    #                 HasLoop2 +=1
    #             else:
    #                 pass

    return
    
    plt.show()

            

                        
           
              

           
              
              


