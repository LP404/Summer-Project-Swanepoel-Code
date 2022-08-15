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

def FancyInterpolate(xInterp,xArray,yArray,AntiNodeX,AntiNodeY,BoolMinima):
    
    #Loc0 = np.where((np.diff(yArray) / max(np.diff(yArray))) == 1)[0][0]
    
    if BoolMinima == True:
        Loc = max(np.where(yArray < AntiNodeY[0])[0])
    else: 
        Loc = np.where(xArray ==  AntiNodeX[0])[0][0]
        
    #Loc = max(np.where((np.around(yArray,2) == 0.50) | (np.around(yArray,2) == 0.49) | (np.around(yArray,2) == 0.51))[0])
    
    #0.5 is currently arbiraty, but frindges don't usually appear until after 50% transmission
    
    x1 = np.append(xArray[0:Loc],AntiNodeX)
    y1 = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(x1)
    x1 = x1[Sort[::1]]
    y1 = y1[Sort[::1]]
        
    yInterp = np.interp(xInterp,x1,y1)
    
    return yInterp


def FindAntiNode(xArray,yArray):
    
    xMax = xArray[argrelextrema(yArray, np.greater)[0]]
    xMin = xArray[argrelextrema(yArray, np.less)[0]]
    
    yMax = yArray[argrelextrema(yArray, np.greater)[0]]
    yMin = yArray[argrelextrema(yArray, np.less)[0]]
    
    
    #Trims front of the maximias
    
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
    
    return xMax, yMax, xMin, yMin

def BackNodeFix(Val1,Val2,xArray,yArray,xMax,yMax,xMin,yMin):
   
   #!!! Make this smart and apply filter when it finds an insonistancy or a very small gap 
   #!!! Currently relies on a magic number of 0.075
   
   
    y = NoiseFilter(Val1,Val2,yArray)
    
    DelMaxLoc = np.where(xMax >= 600)[0]
    DelMinLoc = np.where(xMin >= 600)[0]
    
    yMax = np.delete(yMax,DelMaxLoc)
    xMax = np.delete(xMax,DelMaxLoc)   
    yMin = np.delete(yMin,DelMinLoc)
    xMin = np.delete(xMin,DelMinLoc)
    
    
    xMax1,yMax1,xMin1,yMin1 = FindAntiNode(xArray,y)
    
    
    yMax = np.append(yMax,yMax1[np.where(xMax1 >= 600)[0]])
    xMax = np.append(xMax,xMax1[np.where(xMax1 >= 600)[0]])   
    yMin = np.append(yMin,yMin1[np.where(xMin1 >= 600)[0]])
    xMin = np.append(xMin,xMin1[np.where(xMin1 >= 600)[0]])
    
    return xMax,yMax,xMin,yMin


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

#Note value must be to the nearest 0.5, investigate further
def round_off_rating(Array):
    value = round(Array * 2) / 2
    if value % 1 == 0.5:    
        return value
    
    elif Array % 1 == 0:
         return value + 0.5
        
    else:   
        newValue = value - Array
        if newValue / abs(newValue) == 1:
            return (value - 1) + 0.5
        else:
            return value + 0.5
        


def dFinder(n1,n2,lam1,lam2,M):
    
    if lam1 > lam2:   
        d = (M*lam1*lam2) / (2*((lam1*n2)-(lam2*n1)))
    else:
        d = (M*lam1*lam2) / (2*((lam2*n1)-(lam1*n2)))
    
    
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

def gauss(mean,std,array,yarray):
   for x in range(0,len(array)):
        yarray[x] = (1 / (std * np.sqrt(2*np.pi))) * np.exp(-((array[x]-mean)**2)/(2*std**2))
   return yarray

def main():
    print('Do not run this directly')

    return

if __name__ == "__main__":
    main()
    