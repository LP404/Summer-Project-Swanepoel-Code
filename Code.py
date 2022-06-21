import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#dt = 1

# def FindAntiNode2(xArray,yArray):
#     MaxLocsX = np.array([])
#     MaxLocsY = np.array([])
#     MinLocsX = np.array([])
#     MinLocsY = np.array([])
    
#     for i in range(1,len(xArray)-1):
#         m2ndDer = (yArray[i+1] - (2 * yArray[i]) + yArray[i-1]) / (xArray[i+1] - (2 * xArray[i]) + xArray[i-1])
        
#         if m2ndDer <= 0.1:
#             m1 = (yArray[i] - yArray[i-1]) / (xArray[i] - xArray[i-1])
#             m2 = (yArray[i+1] - yArray[i]) / (xArray[i+1] - xArray[i])
            
#             Trial1 = m1 / abs(m1)
#             Trial2 = m2 / abs(m2)
                            
#             if Trial1 < Trial2:
#                     MaxLocsX = np.append(MaxLocsX,xArray[i])
#                     MaxLocsY = np.append(MaxLocsY,yArray[i])
#             elif Trial2 > Trial1:
#                     MinLocsX = np.append(MinLocsX,xArray[i])
#                     MinLocsY = np.append(MinLocsY,yArray[i])
#             else:
#                 pass
                
                
            


#     return MaxLocsX, MaxLocsY, MinLocsX, MinLocsY

def Trim(xArray,yArray,start,end):
    delPoints = np.where((xArray < start) | (xArray > end))[0]
    
    xArray = np.delete(xArray,delPoints)
    yArray = np.delete(yArray,delPoints)
    
    
    return xArray, yArray



def FancyInterpolate(xInterp,xArray,yArray,AntiNodeX,AntiNodeY):
    
    Loc0 = np.where((np.diff(yArray) / max(np.diff(yArray))) == 1)[0][0]
    Wave0 = xArray[Loc0]
    
    
    Loc = max(np.where(np.around(yArray,0) == 50)[0])
    
    #50 is currently arbiraty, but frindges don't usually appear until after 50% transmission
    
    x1 = np.append(xArray[0:Loc],AntiNodeX)
    y1 = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(x1)
    x1 = x1[Sort[::1]]
    y1 = y1[Sort[::1]]
        
    yInterp = np.interp(xInterp,x1,y1)
    
    return yInterp,Wave0


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

#!!!Once case statments are released, change this
def AntiNodeHighlander(blocks,xArray,yArray,xFull,yFull):
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
        
        
        try:
            FindLoc = np.where(yRevAntiNode < 20)[0][0]
            xRevAntiNode = np.delte(xRevAntiNode,FindLoc)
            yRevAntiNode = np.delte(yRevAntiNode,FindLoc)
        except:
            pass
        
        

    return xRevAntiNode, yRevAntiNode

def FindAntiNode(xArray,yArray,RangeMax):
    
    #RangeMaxNumbers
    #2 is the lowest number
    # Nearest Border - 2 is theoretically the highest but use 15
    #Current setting is arbitary
    
    #!!!Add try catch later
    
    if RangeMax < 2:
        RangeMax = 2
    elif RangeMax > (int(len(xArray)/2) - 2):
        RangeMax = (int(len(xArray)/2) - 2)
    else:
        pass
    
    MaxLocsX = np.array([])
    MaxLocsY = np.array([])
    MinLocsX = np.array([])
    MinLocsY = np.array([])
    
    MinAllowed = False
    MaxAllowed = False
    
    
    #0 is before a potential anti note, 1 is what is being tested to be the antinode, and 2 is the point after the potential antinode
    
    #max(np.where(vars()[files[0]][1] < 1)[0]) Add an engagment mechanisim past 250nm
    #WTF is engement?
    
    
    for i in range(5,len(xArray)-5):
        if i != 0:
            m1 = (yArray[i] - yArray[i-1]) / (xArray[i] - xArray[i-1])
            m2 = (yArray[i+1] - yArray[i]) / (xArray[i+1] - xArray[i])
            

            Trial1 = m1 / abs(m1)
            Trial2 = m2 / abs(m2)
                            
            if Trial1 != Trial2:
                Prev = 0
                Next = 0
                for k in range(1,RangeMax):
                    if i < (len(xArray) - (RangeMax+1)):
                        Prev += (yArray[i-(k-1)] - yArray[i-k]) / (xArray[i-(k-1)] - xArray[i-k])
                        Next += (yArray[i+(k+1)] - yArray[i+k]) / (xArray[i+(k+1)] - xArray[i+k])
                    else:
                        pass
            
                if Prev < 0 and Next > 0:
                    MinAllowed = True
                elif Prev > 0 and Next < 0:
                    MaxAllowed = True
                else:
                    pass
            
                if MaxAllowed == True:
                    MaxLocsX = np.append(MaxLocsX,xArray[i])
                    MaxLocsY = np.append(MaxLocsY,yArray[i])
                    MaxAllowed = False
                elif MinAllowed == True:
                    MinLocsX = np.append(MinLocsX,xArray[i])
                    MinLocsY = np.append(MinLocsY,yArray[i])
                    MinAllowed = False
                else:
                    pass
            else:
                pass

            
    return MaxLocsX, MaxLocsY, MinLocsX, MinLocsY


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
    vars()[files[i]+'T'] = Trim(vars()[files[i]][0],vars()[files[i]][1],xMin,xMax)

#T stands for Truncated
    
    
for i in range(len(files)):
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)], vars()['MinY'+str(i)] = FindAntiNode(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],10)

for i in range(len(files)):
    vars()['MaxBlocks'+str(i)] = BlockFinder(vars()['MaxX'+str(i)])
    vars()['MinBlocks'+str(i)] = BlockFinder(vars()['MinX'+str(i)])

for i in range(len(files)):
    vars()['xNewMax'+str(i)], vars()['yNewMax'+str(i)] = AntiNodeHighlander(vars()['MaxBlocks'+str(i)],vars()['MaxX'+str(i)], vars()['MaxY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])
    vars()['xNewMin'+str(i)], vars()['yNewMin'+str(i)] = AntiNodeHighlander(vars()['MinBlocks'+str(i)],vars()['MinX'+str(i)], vars()['MinY'+str(i)],vars()[files[i]][0],vars()[files[i]][1])

#All the xValues are the same, only need to do this once
xP = np.linspace(min(vars()[files[0]+'T'][0]),max(vars()[files[0]+'T'][0]),10001)

for i in range(len(files)):
    vars()['yPMax'+str(i)], vars()['MaxLocWave'+str(i)]  = FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMax'+str(i)],vars()['yNewMax'+str(i)])
    vars()['yPMin'+str(i)], vars()['MinLocWave'+str(i)] = FancyInterpolate(xP,vars()[files[i]+'T'][0],vars()[files[i]+'T'][1],vars()['xNewMin'+str(i)],vars()['yNewMin'+str(i)])






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
    
    
