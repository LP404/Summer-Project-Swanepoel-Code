import numpy as np
import scipy as sp





# def FindAntiNode2(xArray,yArray,dt):
#     MaxLocsX = np.array([])
#     MaxLocsY = np.array([])
#     MinLocsX = np.array([])
#     MinLocsY = np.array([])
    
#     for i in range(1,len(xArray)-1):
#         m2ndDer = (yArray[i+1] - (2 * yArray[i]) + yArray[i-1]) / dt**2
        
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

def Trim(Array,start,end):
    delPoints = np.where((Array[0] < start) | (Array[0] > end))[0]
    
    
    newArray = np.array([np.delete(Array[0],delPoints),np.delete(Array[1],delPoints)])
    
    return newArray



def FancyInterpolate(xInterp,xArray,yArray,AntiNodeX,AntiNodeY):
    
    #Loc0 = np.where((np.diff(yArray) / max(np.diff(yArray))) == 1)[0][0]
    
    Loc = max(np.where(np.around(yArray,2) == 0.50)[0])
    
    #0.5 is currently arbiraty, but frindges don't usually appear until after 50% transmission
    
    x1 = np.append(xArray[0:Loc],AntiNodeX)
    y1 = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(x1)
    x1 = x1[Sort[::1]]
    y1 = y1[Sort[::1]]
        
    yInterp = np.interp(xInterp,x1,y1)
    
    return yInterp


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


def ArrayInterlace(Array1, Array2):
    Combi = np.empty(Array1.size + Array2.size, dtype=Array1.dtype)
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
    
    if min(xMax) <= min(xMin):
        Combi = ArrayInterlace(xMax, (-1 * xMin))     
    else:
        Combi = ArrayInterlace((-1* xMin), * xMax)
        
    for i in range(1,len(Combi)):
        if (Combi[i-1]/Combi[i-1]) != (Combi[i]/Combi[i]):
            pass
        elif (Combi[i-1]/abs(Combi[i-1])) == (Combi[i]/abs(Combi[i])) == 1:  
            Loc = np.where(Combi[i-1] == xMax)[0][0]
            Loc1 = np.where(Combi[i] == xMax)[0][0]
            
            if yMax[Loc] <= yMax[Loc1]:
                yMax = np.delete(yMax,Loc)
                xMax = np.delete(yMax,Loc)
            else:
                yMax = np.delete(yMax,Loc1)
                xMax = np.delete(yMax,Loc1)   
                  
        elif (Combi[i-1]/abs(Combi[i-1])) == (Combi[i]/abs(Combi[i])) == -1:
            Loc = np.where(abs(Combi[i-1]) == xMin)[0][0]
            Loc1 = np.where(abs(Combi[i]) == xMin)[0][0]
            if yMin[Loc] >= yMin[Loc1]:                    
                yMin = np.delete(yMin,Loc)
                xMin = np.delete(yMin,Loc)
            else:
                yMin = np.delete(yMin,Loc1)
                xMin = np.delete(yMin ,Loc1) 
        else:
            pass
                
                    
    return xMax, yMax, xMin, yMin


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
    
    #These need to be at 5 and I don't remember why
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

def DYThorLabs(Array):
    
      newXarray = np.array([])
      newYarray = np.array([])
    
      newXarray = np.append(newXarray,Array[0][0])
      backTotalCount = 0
      backTotal = 0
    
      for i in range(1,len(Array[0])):
          if Array[0][i] == Array[0][i-1]:
             backTotalCount += 1
          else:   
             newXarray = np.append(newXarray,Array[0][i])
             newYarray = np.append(newYarray,np.mean(Array[1][backTotal:i]))
             backTotal = backTotalCount + 1
    
    
      backTotal = backTotalCount + 1
      newYarray = np.append(newYarray,np.mean(Array[1][backTotal:i]))
    
      FixedArray = np.array([newXarray,newYarray])
        
      return FixedArray

def SecondTrim(xMax,yMax,xMin,yMin):

    if xMax[-1] > xMin[-1]:
        delPointsMax = np.array([0,1,len(xMax)-2,len(xMax)-1])
        delPointsMin = np.array([0,1,len(xMin)-1])
    else:
        delPointsMax = np.array([0,1,len(xMax)-2,len(xMax)-1])
        delPointsMin = np.array([0,1,len(xMin)-2,len(xMin)-1])        
        
    
    xMaxT = np.delete(xMax,delPointsMax)
    yMaxT = np.delete(yMax,delPointsMax)
    xMinT = np.delete(xMin,delPointsMin)
    yMinT = np.delete(yMin,delPointsMin)
    

    return xMaxT, yMaxT, xMinT, yMinT


def FindNearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
    
    
def nFinder(TM,Tm,lamb,s):
    
    N = (2*s) * ((TM - Tm) / (TM*Tm)) + ((s**2 + 1)/2)
    
    nA = np.sqrt((N**2 - s**2))
    n = np.sqrt((N + nA)) 
    
    
    return n

def round_off_rating(Array):
    return round(Array * 2) / 2


def dFinder(n1,n2,lam1,lam2):
    
    if lam1 > lam2:   
        d = (lam1*lam2) / (2*((lam1*n2)-(lam2*n1)))
    else:
        d = (lam1*lam2) / (2*((lam2*n1)-(lam1*n2)))
    
    
    return d

def xFinder(n,s,TM,Tm):
    
    F = 4* n**2 * s * ((TM+Tm)/(TM*Tm))    
    
    Num = F - np.sqrt( F**2 - ((n**2 - 1)**3 * (n**2 - s**4))   )
    Dem = (n-1)**3 * (n - s**2)
    
    x = Num/Dem
    
    return x

def Sub_nFinder(T):
    ns = 1/T + np.sqrt(((1/T)**2 - 1))
    
    return ns

def main():
    print('Do not run this directly')

    return

if __name__ == "__main__":
    main()
    