import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#dt = 1

def Trim(xArray,yArray,start,end):
    delPoints = np.where((xArray < start) | (xArray > end))[0]
    
    xArray = np.delete(xArray,delPoints)
    yArray = np.delete(yArray,delPoints)
    
    
    return xArray, yArray


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



def DistanceChecker(Array):
    BlockLengths = np.array([])
    PreviousGap = Array[1] - Array[0]
    BlockCount = 0
    for i in range(2,len(Array)):
        Gap = abs(Array[i] - Array[i-1])
        if Gap >= 3*PreviousGap :
            BlockLengths = np.append(BlockLengths,BlockCount)
            BlockCount = 0
        else:
            BlockCount += 1
            
            PreviousGap = Gap
            Gap = 0
    
    
    return BlockLengths

def FindAntiNode(xArray,yArray):
    MaxLocsX = np.array([])
    MaxLocsY = np.array([])
    MinLocsX = np.array([])
    MinLocsY = np.array([])
    
    MinAllowed = False
    MaxAllowed = False
    #2 is the lowest number
    RangeMax = 2
    
    #0 is before a potential anti note, 1 is what is being tested to be the antinode, and 2 is the point after the potential antinode
    
    #max(np.where(vars()[files[0]][1] < 1)[0]) Add an engagment mechanisim past 250nm
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



path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\Data'))


#The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
for i in range(len(files)):
    files[i] = files[i].rstrip(".txt")

#The next three loops import and assign varable names to the data 
for i in range(len(files)):
    vars()[files[i]] = np.loadtxt(open(path + "\\" + files[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T

for i in range(len(files)):
    vars()[files[i]+'T'] = Trim(vars()[files[i]][0],vars()[files[i]][1],200,850)

#T stands for Truncated
    
    
for i in range(len(files)):
    vars()['MaxX'+str(i)], vars()['MaxY'+str(i)], vars()['MinX'+str(i)], vars()['MinY'+str(i)] = FindAntiNode(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1])

#All the xValues are the same, only need to do this once
xP = np.linspace(min(vars()[files[0]+'T'][0]),max(vars()[files[0]+'T'][0]),10001)

for i in range(len(files)):
    
    vars()['yPMaxs'+str(i)] = np.interp(xP,vars()['MaxX'+str(i)],vars()['MaxY'+str(i)])
    vars()['yPMins'+str(i)] = np.interp(xP,vars()['MinX'+str(i)],vars()['MinY'+str(i)])
    
    



for i in range(len(files)):
    plt.figure(i)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='darkgray', linestyle='--')
    plt.plot(vars()[files[i]+'T'][0],vars()[files[i]+'T'][1])
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Transmission/Absorbance")
    # plt.scatter(vars()['MaxX'+str(i)],vars()['MaxY'+str(i)], color = 'black', marker = "x")
    # plt.scatter(vars()['MinX'+str(i)],vars()['MinY'+str(i)], color = 'red', marker = "x")
    # plt.plot(xP,vars()['yPMaxs'+str(i)], color = 'black', linestyle="dashed")
    # plt.plot(xP,vars()['yPMins'+str(i)], color = 'red', linestyle="dashed")
    
    
