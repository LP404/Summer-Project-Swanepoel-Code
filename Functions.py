import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter
from scipy import stats
from scipy.optimize import curve_fit

def ArrayCondenser(x,y):
      newXarray = np.array([])
      newYarray = np.array([])
    
      newXarray = np.append(newXarray,x[0])
      backTotalCount = 0
      backTotal = 0
    
      for i in range(1,len(x)):
          if x[i] == x[i-1]:
             backTotalCount += 1
          else:   
             newXarray = np.append(newXarray,x[i])
             newYarray = np.append(newYarray,np.mean(y[backTotal:i]))
             backTotal = backTotalCount + 1
    
    
      backTotal = backTotalCount + 1
      newYarray = np.append(newYarray,np.mean(y[backTotal:i]))
    
      return newXarray,newYarray

def ArrayFix(x,y):
    
    x2 = np.array([x[0]])
    y2 = np.array([y[0]])
    
    for i in range(1,len(x)):
        if x[i] != x[i-1]:
            x2 = np.append(x2,x[i])
            y2 = np.append(y2,y[i])
        else:
            yHold = (y[i] + y[i-1]) / 2
            y2[i-1] = yHold
    
    return x2, y2

def LineValFinder(xArray,yArray,guage,setting,xMin,xMax,inverse):
    
    if xMin > max(xArray) or xMax < min(xArray):
        return print('Error incorrect Max/Min values max: '+str(max(xArray))+' vs stated of '+str(xMax)+' and min: '+str(min(xArray))+' vs stated of '+str(xMin))
    else:   
        if setting == 'min':
            point = np.where(np.gradient(yArray) == min(np.gradient(yArray)))[0][0]      
        elif setting == 'bg' and inverse == False:
            allowed = np.where((xArray < xMax) & (xArray > xMin))[0]
            xNew = xArray[allowed]
            yNew = yArray[allowed]
            point = np.where(np.gradient(yNew) == min(np.gradient(yNew)))[0][0]  
        elif setting == 'bg' and inverse == True:
            allowed = np.where((xArray < xMax) & (xArray > xMin))[0]
            xNew = xArray[allowed]
            yNew = yArray[allowed]
            #Due to the conversion to eV, the graph is effectivley "reversed, i.e +/- graidents switch in value)
            point = np.where(-np.gradient(yNew) == min(-np.gradient(yNew)))[0][0]  
        else:
            point = np.where(np.gradient(yArray) == max(np.gradient(yArray)))[0][0]
            
        xSec = xNew[point - guage : point + guage]
        ySec = yNew[point - guage : point + guage] 
         
         
        m,c,r,p,err = stats.linregress(xSec,ySec)
        
        
         
        return m,c,err

def LineIntersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


#Takes in two arrayes and finds all indices where one is in another
#Assumes each value in the bigArray is unique
def ValueFinder(guageArray,bigArray):
    num = len(guageArray)
    finalArrayInd = np.zeros(num)
    for i in range(len(guageArray)):
        guess = guageArray[i]
        point = np.where(guess == bigArray)[0][0]
        finalArrayInd[i] = int(point)
        
        
    finalArrayInd = finalArrayInd.astype(int)
    return finalArrayInd


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

def Richards(x,A,K,B,v,Q,C):
    num = K - A
    denA = C + (Q*(np.exp((-B*x))))
    denB = denA**(1/v)
    frac = num/denB
    final = A + frac
    
    return final

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

#For use with the ThorLabsSapphire.txt
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

def FindNearestVal(array, value):
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

def gauss(mean,std,array):
   yarray = np.zeros(len(array)) 
   for x in range(0,len(array)):
        yarray[x] = (1 / (std * np.sqrt(2*np.pi))) * np.exp(-((array[x]-mean)**2)/(2*std**2))
   return yarray

def ThicknessAcceptance(Lambda,d,CutOff,dErr,stdScale):

    LambdaInds = Lambda.argsort()
    d = d[LambdaInds[::1]]
    Lambda = Lambda[LambdaInds[::1]]
    dErr = dErr[LambdaInds[::1]]

    LambdaDiscard = Lambda[0:CutOff]
    
    d = d[CutOff:]
    Lambda = Lambda[CutOff:]
    dErr = dErr[CutOff:]

    dInds = d.argsort()
    d = d[dInds[::1]]
    Lambda = Lambda[dInds[::1]]
    dErr = dErr[dInds[::1]]
    
    
    dStd = np.std(d)
    dMean = np.mean(d)
    
    dGauss = gauss(dMean,stdScale*dStd,d)
    
    Multiplier = max(dGauss) / max(Lambda)
    
    GoodInds = np.where(Multiplier*Lambda < dGauss)
    BadInds = np.where(Multiplier*Lambda >= dGauss)
    
                
    RejectedLambda = Lambda[BadInds]
    RejectedLambda = np.append(RejectedLambda,LambdaDiscard)
    Newd = d[GoodInds]
    
    NewdErr = dErr[GoodInds]
    dErr = avgUncert(NewdErr)
    
    newMean = np.mean(Newd)
    error = np.std(Newd) / np.sqrt(len(Newd))
    
    
    return newMean,error,RejectedLambda,dErr


def PerToAbs(PerUncert,Val):
    AbsUncert = (PerUncert/100)*Val
    return AbsUncert

def AbsToPer(AbsUncert,Val):
    PerUncert = (AbsUncert/Val)*100
    return PerUncert

#Used to add uncetanties together for either absolutes for add/subtract or percentages for multuplt/deivide
def LinearCombine(Uncert1,Uncert2):
    CombiUncert = np.sqrt((abs(Uncert1)**2 + abs(Uncert2)**2))
    return CombiUncert

def MultiDivUncert(Uncert1,Uncert2,Val1,Val2):
    CombiUncert = (np.sqrt(((abs(Uncert1)/Val1)**2) + ((abs(Uncert2)/Val2)**2))) * 100
    return CombiUncert

def ExpoUncert(Uncert,Val,Expo):
    NewUncert = Expo*((abs(Uncert)/Val)*100)
    return NewUncert

def ScalingUncert(Uncert,Scalar):
    return (Uncert*Scalar)

def nUncertFinder(TUncert,TM,Tm,n):
    
    NUncertNume = LinearCombine(TUncert,TUncert)
    NUncertDenom = MultiDivUncert(TUncert,TUncert,TM,Tm)
    
    Tdiff = TM-Tm
    Tprod = TM*Tm
    
    NUncertNumePer = AbsToPer(NUncertNume,Tdiff)
    
    NUncertPer = LinearCombine(NUncertNumePer,NUncertDenom)
    
    NComp = Tdiff/Tprod
    
    NUncertAbs = PerToAbs(NUncertPer,NComp)
    
    NUncert = ScalingUncert(NUncertAbs,2)
    
    nUncertPer = ExpoUncert(NUncert,n,0.5)
    
    nUncert = PerToAbs(nUncertPer,n)

    return nUncert


def dUncert(Lam1,Lam2,n1,n2,d,LamUncert,n1Uncert,n2Uncert):
    
    LamProd = Lam1*Lam2
    
    Prod12 = Lam1*n2
    Prod21 = Lam2*n1
    
    
    dDom = 2*abs((Lam1*n2 - Lam2*n1))
    
    dNumUncert = MultiDivUncert(LamUncert, LamUncert, Lam1, Lam2)
    
    dNumUncertAbs = PerToAbs(dNumUncert,LamProd)
    
    dDomUncertPartA = MultiDivUncert(LamUncert,n1Uncert,Lam2,n1) 
    dDomUncertPartB = MultiDivUncert(LamUncert,n2Uncert,Lam1,n2)
    
    dDomUncertPartAabs = PerToAbs(dDomUncertPartA,Prod21)
    dDomUncertPartBabs = PerToAbs(dDomUncertPartB,Prod12)
    
    dDomUncertPartC = LinearCombine(dDomUncertPartAabs,dDomUncertPartBabs)
    
    dDomUncertPartD = ScalingUncert(dDomUncertPartC, 2)
    
    dUncert = MultiDivUncert(dNumUncertAbs,dDomUncertPartD,LamProd,dDom)
    
    dUncertAbs = PerToAbs(dUncert, d)
    
    return dUncertAbs


def dUncerIgnN(Lam1,Lam2,n1,n2,d,LamUncert):
    
    LamProd = Lam1*Lam2
    dNumUncert = MultiDivUncert(LamUncert, LamUncert, Lam1, Lam2)
    dNumUncertAbs = PerToAbs(dNumUncert,LamProd)
    
    dDom = 2*abs((Lam1*n2 - Lam2*n1))
    
    dDomUncertPartA = LinearCombine(LamUncert,LamUncert)
    dDomUncertPartB = ScalingUncert(dDomUncertPartA, 2)
    dUncert = MultiDivUncert(dNumUncertAbs,dDomUncertPartB,LamProd,dDom)
    
    dUncertAbs = PerToAbs(dUncert, d)
    
    
    return dUncertAbs
    

def avgUncert(vals):
    
    Num = 0
    
    for i in vals:
        Num += i**2
        
    Num = np.sqrt(Num)
    
    Den = len(vals)
    
    Uncert = Num/Den
    
    return Uncert


def mUncert(n,d,Lam,nUncert,dUncert,LamUncert,m):
    
    Num = n*d
    
    NumUncert = MultiDivUncert(nUncert,dUncert,n,d)
    
    NumUncertAbs = PerToAbs(NumUncert,Num) 
    
    mUncertPer = MultiDivUncert(NumUncertAbs,LamUncert,Num,Lam)
    
    mUncertAbs = PerToAbs(mUncertPer,m)
    
    mUncert = 2*mUncertAbs
    
    return mUncert


def d2Uncert(n,m,Lam,nUncert,mUncert,LamUncert,d):

    Num = m*Lam
    
    NumUncert = MultiDivUncert(mUncert,LamUncert,m,Lam)
    
    NumUncertAbs = PerToAbs(NumUncert,Num) 
    
    dUncertPer = MultiDivUncert(NumUncertAbs,nUncert,Num,n)
    
    dUncertAbs = PerToAbs(dUncertPer,d)
    
    dUncert = dUncertAbs/2    
    
    return dUncert


def n2Uncert(d,m,Lam,dUncert,mUncert,LamUncert,n):
    
    Num = m*Lam
    
    NumUncert = MultiDivUncert(mUncert,LamUncert,m,Lam)
    
    NumUncertAbs = PerToAbs(NumUncert,Num) 
    
    nUncertPer = MultiDivUncert(NumUncertAbs,dUncert,Num,d)
    
    nUncertAbs = PerToAbs(nUncertPer,n)
    
    nUncert = nUncertAbs/2    
    
    return nUncert

def ellipseParamFinder(x1,y1,x2,y2):
    line1 = [[200,y1] , [800,y1]]
    line2 = [[x2,0] , [x2,1]]
    x0, y0 = LineIntersection(line1 ,line2)
    a = abs(x1-x0)
    b = abs(y2-y0)
    c = np.sqrt((a**2 + b**2))
    
    return x0,y0,a,b,c

def quarterEllipseFinder(x,a,b,x0,y0):
    y = y0 + (((b/a) * (np.sqrt((a**2-((x-x0)**2))))))
    return y
    
def logCurve(x,a,b,c):
    
    y = np.emath.logn(a,x-c) + b
        
    return y 

#Note define an x before passing in
def LineComp(x,m1,m2,c1,c2,delta):
    m1Index = 0
    m2Index = 0
    MaxOccur = 0
    for i in range(len(m1)):
        y1 = m1[i]*x + c1[i]
        for j in range(len(m2)):
            y2 = m2[j]*x + c2[j]
            checkArray = y1/y2
            for k in range(len(checkArray)):
                if checkArray[k] < 1:
                    checkArray[k] = (1 / checkArray[k])
                else:
                    pass
                
            oneOccur = np.count_nonzero(np.around(checkArray,delta) == 1)
            if oneOccur > MaxOccur:
                m1Index = i
                m2Index = j
                MaxOccur = oneOccur
            else:
                pass
    return MaxOccur, m1Index, m2Index

def Reinterp(xInterp,xArray,yArray,AntiNodeX,AntiNodeY,CutOff):
    #Since we know the cutoff is 0.6
    #Given we know how the array is constructed we can assume the next index is for the first maxima/minima

    Loc = max(np.where((np.around(yArray,2) == CutOff) | (np.around(yArray,2) == CutOff-0.01) | (np.around(yArray,2) == CutOff+0.01))[0])
    
    
    xNuvo = np.append(xArray[0:Loc],AntiNodeX)
    yNuvo = np.append(yArray[0:Loc],AntiNodeY)
        
    Sort = np.argsort(xNuvo)
    xNuvo = xNuvo[Sort[::1]]
    yNuvo = yNuvo[Sort[::1]]


    yNuvoInd = np.where(yNuvo == FindNearestVal(yNuvo,CutOff))[0][0]
    y1 = yNuvo[yNuvoInd - 4]
    x1 = xNuvo[yNuvoInd - 4]
    y2 = yNuvo[yNuvoInd + 1]
    x2 = xNuvo[yNuvoInd + 1]
    x0, y0, a, b, c = ellipseParamFinder(x1,y1,x2,y2)   
    
    xEllipse = np.arange(x1, x2+1, 1)
    yEllipse = quarterEllipseFinder(xEllipse,a,b,x0,y0)
        
    guess = np.array([np.exp(1),1,260])
    

    yTopStart =   yNuvo[yNuvoInd+1:]
    xTopStart =   xNuvo[yNuvoInd+1:]
    
    
    TransFit, Resi = curve_fit(logCurve,xTopStart, yTopStart,guess, maxfev=500000000)
       
    xTopStart2 = np.arange(xNuvo[yNuvoInd +1], 801, 1)
        
    yFitbutIDK = logCurve(xTopStart2 , TransFit[0], TransFit[1], TransFit[2])

    
    
    
    
    #Since step size is 1 we can use np.gradient
    
    #!!! Note: Change X values to cover all x points for final fit
    

    gradEllipse = np.gradient(yEllipse)
    bEllipse = yEllipse - (gradEllipse*xEllipse)
    
    gradyFitbutIDK = np.gradient(yFitbutIDK)
    bFitbutIDK = yFitbutIDK - (gradyFitbutIDK*xTopStart2)
 
    allowance = 20

    xArr = np.linspace(x0-allowance,x0+allowance,100001)
    IncLoc = LineComp(xArr,gradEllipse,gradyFitbutIDK,bEllipse,bFitbutIDK,5)



    yArr = (xArr*gradEllipse[IncLoc[1]] + bEllipse[IncLoc[1]])
    yArr1 = (xArr*gradyFitbutIDK[IncLoc[2]] + bFitbutIDK[IncLoc[2]])
    
    xArr = np.round(xArr).astype(int)

    GeoMeanNewyArr = np.sqrt((yArr * yArr1))
    

    newxArr, newGeoMeanNewyArr = ArrayCondenser(xArr,GeoMeanNewyArr)

    xStitched,yStitched = Stitcher(x0,xEllipse,yEllipse,xTopStart2,yFitbutIDK,newxArr,newGeoMeanNewyArr,xArray,yArray)
    xStitchedFixed,yStitchedFixed = ArrayFix(xStitched,yStitched)
   
    yInterp = np.interp(xInterp,xStitchedFixed,yStitchedFixed)
   
    return yInterp

def Stitcher(x0,xElip,yElip,xCurveFit,yCurveFit,xLineFit,yLineFit,xOriginal,yOriginal):
        
    Cut = np.where(xLineFit == x0)[0][0]
        
    xLine1 = xLineFit[:Cut]
    xLine2 = xLineFit[Cut:]
    
    yLine1 = yLineFit[:Cut]
    yLine2 = yLineFit[Cut:]    
    
    LenVal1 = len(xElip)
    LenVal2 = len(xCurveFit)
    
    ActualVal1 = len(xLine1)
    ActualVal2 = len(xLine2)
    
    PadVal1 = abs(LenVal1 - ActualVal1)
    PadVal2 = abs(LenVal2 - ActualVal2)
    
    if LenVal1 < ActualVal1:
    
        xPaddedElipAlt = np.append(np.zeros(PadVal1),xElip)
        yPaddedElipAlt = np.append(np.zeros(PadVal1),yElip)
    
        MinValElipIndAlt = np.argmin(abs(yLine1 - yPaddedElipAlt))
    
        NewXElipAlt = xLine1[MinValElipIndAlt:]
        NewYElipAlt = yLine1[MinValElipIndAlt:]
        
        NewXElip2Alt = xPaddedElipAlt[:MinValElipIndAlt]
        NewYElip2Alt = yPaddedElipAlt[:MinValElipIndAlt]
    
        NewXElip2Alt = np.trim_zeros(NewXElip2Alt)
        NewYElip2Alt = np.trim_zeros(NewYElip2Alt) 
        
        FinalXInterA = np.append(NewXElip2Alt,NewXElipAlt)
        FinalYInterA = np.append(NewYElip2Alt,NewYElipAlt)
    
    else:
        xPaddedElip = np.append(np.zeros(PadVal1),xLine1)
        yPaddedElip = np.append(np.zeros(PadVal1),yLine1)
    
        MinValElipInd = np.argmin(abs(yElip-yPaddedElip))
    
        NewXElip = xPaddedElip[MinValElipInd:]
        NewYElip = yPaddedElip[MinValElipInd:]
        
        NewXElip2 = xElip[:MinValElipInd]
        NewYElip2 = yElip[:MinValElipInd]    
    
        FinalXInterA = np.append(NewXElip2,NewXElip)
        FinalYInterA = np.append(NewYElip2,NewYElip)
    
    
    
    if LenVal2 < ActualVal2:
    
        xPaddedCurveAlt = np.append(xCurveFit,np.zeros(PadVal2))
        yPaddedCurveAlt = np.append(yCurveFit,np.zeros(PadVal2))
            
        
        MinValCurveInd = np.argmin(abs(yLine2-yPaddedCurveAlt))
        
        NewXCurveAlt = xLine2[:MinValCurveInd]
        NewYCurveAlt = yLine2[:MinValCurveInd]  
        
        NewXCurve2Alt = xPaddedCurveAlt[MinValCurveInd:]
        NewYCurve2Alt = yPaddedCurveAlt[MinValCurveInd:]      
        
        NewXCurve2Alt = np.trim_zeros(NewXCurve2Alt)
        NewYCurve2Alt = np.trim_zeros(NewYCurve2Alt)   
        
        FinalXInterB = np.append(FinalXInterA,NewXCurveAlt)
        FinalXInterC = np.append(FinalXInterB,NewXCurve2Alt)
    
        FinalYInterB = np.append(FinalYInterA,NewYCurveAlt)
        FinalYInterC = np.append(FinalYInterB,NewYCurve2Alt)
        
    else: 
    
        xPaddedCurve = np.append(xLine2,np.zeros(PadVal2))
        yPaddedCurve = np.append(yLine2,np.zeros(PadVal2))
            
        
        MinValCurveInd = np.argmin(abs(yCurveFit-yPaddedCurve))
        
        NewXCurve = xPaddedCurve[:MinValCurveInd]
        NewYCurve = yPaddedCurve[:MinValCurveInd]  
        
        NewXCurve2 = xCurveFit[MinValCurveInd:]
        NewYCurve2 = yCurveFit[MinValCurveInd:]   
        
        FinalXInterB = np.append(FinalXInterA,NewXCurve)
        FinalXInterC = np.append(FinalXInterB,NewXCurve2)
    
        FinalYInterB = np.append(FinalYInterA,NewYCurve)
        FinalYInterC = np.append(FinalYInterB,NewYCurve2)
    
    
    startInds = np.where(xOriginal < min(FinalXInterC))[0]
    FinalX = np.append(xOriginal[startInds],FinalXInterC)
    FinalY = np.append(yOriginal[startInds],FinalYInterC)


    return FinalX, FinalY

def main():
    print('Do not run this directly')

    return

if __name__ == "__main__":
    main()
    