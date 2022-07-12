import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter




def Trim(Array,start,end):
    delPoints = np.where((Array[0] < start) | (Array[0] > end))[0]
    
    
    newArray = np.array([np.delete(Array[0],delPoints),np.delete(Array[1],delPoints)])
    
    return newArray

def NoiseFilter(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    zi = lfilter_zi(b, a)
    z, _ = lfilter(b, a, yArray, zi=zi*yArray[0])
    z2, _ = lfilter(b, a, z, zi=zi*z[0])
    yFiltered = filtfilt(b, a, yArray)
    
    return yFiltered

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


def FindAntiNode(xArray,yArray):
    
    Xmaxima = xArray[argrelextrema(yArray, np.greater)[0]]
    Xminima = xArray[argrelextrema(yArray, np.less)[0]]
    
    Ymaxima = yArray[argrelextrema(yArray, np.greater)[0]]
    Yminima = yArray[argrelextrema(yArray, np.less)[0]]
    
    return Xmaxima, Ymaxima, Xminima, Yminima

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
    