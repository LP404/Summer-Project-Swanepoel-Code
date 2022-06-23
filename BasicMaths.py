import numpy as np

#Factor in substrate
#Add in reflection?

def nFinder(TM,Tm,lamb,s):
    
    N = (2*s) * ((TM - Tm) / (TM*Tm)) + ((s**2 + 1)/2)
    
    nA = np.sqrt((N**2 - s**2))
    n = np.sqrt((N + nA)) 
    
    
    return n

def dFinder(n1,n2,lam1,lam2):
    
    d = (lam1*lam2) / (2*((lam1*n2)-(lam2*n1)))
    
    return d


def xFinder(n,s,TM,Tm):
    
    F = 4* n**2 * s * ((TM+Tm)/(TM*Tm))    
    
    Num = F - np.sqrt( F**2 - ((n**2 - 1)**3 * (n**2 - s**4))   )
    Dem = (n-1)**3 * (n - s**2)
    
    x = Num/Dem
    
    return x


#Maxs

n1 = nFinder(62.3486,59.6286,446,1.75)
n2 = nFinder(63.6930,60.8116,534,1.75)

d = dFinder(n2,n1,534,446)

x1 = xFinder(n1,1.78,62.3486,59.6286)
x2 = xFinder(n2,1.775,63.6930,60.8116)

k1 = (-446 / (4 * np.pi * d)) * np.log(x1)
k2 = (-534 / (4 * np.pi * d)) * np.log(x2)

#Mins
n3 = nFinder(62.8434,60.1293,480,1.76)
n4 = nFinder(64.2820,61.4183,582,1.755)

d1 = dFinder(n4,n3,582,480)

x3 = xFinder(n3,1.76,62.8434,60.1293)
x4 = xFinder(n4,1.765,64.2820,61.4183)

k3 = (-480 / (4 * np.pi * d1)) * np.log(x3)
k4 = (-582 / (4 * np.pi * d1)) * np.log(x4)
