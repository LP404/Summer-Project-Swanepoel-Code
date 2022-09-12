import numpy as np

def dFinder(n1,n2,lam1,lam2,M):
    
    if lam1 > lam2:   
        d = (M*lam1*lam2) / (2*((lam1*n2)-(lam2*n1)))
    else:
        d = (M*lam1*lam2) / (2*((lam2*n1)-(lam1*n2)))
    
    
    return d

def nFinder(TM,Tm,s):
    
    N = (2*s) * ((TM - Tm) / (TM*Tm)) + ((s**2 + 1)/2)
    
    nA = np.sqrt((N**2 - s**2))
    n = np.sqrt((N + nA)) 
    
    
    return n


LamList = np.array([286,313,350,402,475,596,786])
TMList = np.array([0.583027823,0.609188692,0.620723579,0.630604666,0.64205607,0.654810598,0.685658897])
TmList = np.array([0.567578925,0.579432267,0.592258056,0.603419348,0.61491902,0.631710221,0.640391554])
sList = np.array([1.9403124596558827,1.8976911029453674,1.866949708602932,1.840751273585143,1.8107335892234522,1.7884435041583395,1.7702315367671904])
dList = np.array([])
dList2 = np.array([])


NList = ((2*sList) * ((TMList-TmList)/(TMList*TmList))) + ((sList**2 + 1)/2)

nList = np.sqrt((NList + np.sqrt((NList**2 - sList**2))))

for i in range(1,len(LamList)):
    dVal2 = dFinder(nList[i],nList[i-1],LamList[i],LamList[i-1],1)
    dList2 = np.append(dList2,dVal2)

dSpecial2 = dFinder(nList[len(LamList)-1],nList[0],LamList[len(LamList)-1],LamList[0],(len(LamList)-1))
dList2 = np.append(dList2,dSpecial2)

mList = np.around(((np.mean(dList2[1:]) * 2 * nList) / LamList),0)
d2List2 = (mList * LamList) / (2 * nList)

n2List2 = (mList * LamList) / (2 * np.mean(d2List2))