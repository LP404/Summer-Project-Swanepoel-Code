import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import os
import csv


def Acalc(n,s):
    return 16 * (n**2) * s

def Bcalc(n,s):
    
    P1 = (n + 1)**3
    P2 = (n + s**2)
    
    return  P1 * P2

def Ccalc(n,s):
    
    P1 = (n - 1)**2
    P2 = ((n**2) + (s**2))
    
    return  2 * P1 * P2
    
def Dcalc(n,s):
    
    P1 = (n - 1)**3
    P2 = (n - s**2)
    
    return  P1 * P2

#A1 and A2 calc is the same

def A1calc(n,s,k):
    return 16 * s * ((n**2) + (k**2))

def B1calc(n,s,k):
    
    P1 = ((n + 1)**2) + k**2
    P2 = ((n + 1) * (n + (s**2))) + k**2
    
    return  P1 * P2

def C1calc(n,s,k,phi):
    
    P1 = ((n**2) - 1 + (k**2))
    P2 = ((n**2) - (s**2) + (k**2))
    P3 = 2*(k**2) * ((s**2) + 1)
    
    P4 = ((P1*P2) - P3) * 2*np.cos(phi)
    
    P5 = 2 * ((n**2) - (s**2) + (k**2))
    P6 = ((s**2) + 1)
    P7 = ((n**2) - 1 + (k**2))
    
    P8 = -1 * k * (P5 + (P6*P7)) * 2 * np.sin(phi)
    
    return P4 - P8


def D1calc(n,s,k):
    
    P1 = ((n - 1)**2) + k**2
    P2 = ((n - 1) * (n - (s**2))) + k**2
    
    return  P1 * P2

def B2calc(n,s,k):
    
    P1 = ((n + 1)**2) + k**2
    P2 = ((n + (s**2))) + k**2
    
    return  P1 * P2

def C2calc(n,s,k,phi):
    
    P1 = ((n**2) - 1 + (k**2))
    P2 = ((n**2) - (s**2) + (k**2))
    P3 = 4 * (k**2) * s
    
    P4 = ((P1*P2) - P3) * 2*np.cos(phi)
    
    P5 = 2 * ((n**2) - (s**2) + (k**2))
    P6 = 2 * s * ((n**2) - 1 + (k**2))
    
    P8 = -1 * k * (P5 + P6) * 2 * np.sin(phi)
    
    return P4 - P8

def D2calc(n,s,k):
    
    P1 = ((n - 1)**2) + k**2
    P2 = ((n - (s**2))) + k**2
    
    return  P1 * P2


def phiCalc(n,d,lam):
    
    P1 = n * d
    P2 = P1/lam
    
    return 4*np.pi*P2

def aCalc(k,lam):
    
    P1 = 4*np.pi*k
    
    return P1/lam

def xCalc (a,d):
    return np.exp((-1*a*d))

def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]


def ktoa(k,lam):
    return (4*np.pi * k)/lam



d = 351.7

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\Absorption'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))



suffix = ".csv"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
    rows = []

    with open(str(path)+"\\"+files[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header = next(csvreader)
    
        for row in csvreader:
            rows.append(row)
    
    vars()[f'{files[i]} Data'] = rows
    
    for j in range(len(vars()[f'{files[i]} Data'])):
        for k in range(len(vars()[f'{files[i]} Data'][j])):
            try:
                vars()[f'{files[i]} Data'][j][k] = float(vars()[f'{files[i]} Data'][j][k])
            except ValueError:
                pass

xP = np.array(ListExtract(vars()[f'{files[i]} Data'],0))
Traw = np.array(ListExtract(vars()[f'{files[i]} Data'],1))
Tm = np.array(ListExtract(vars()[f'{files[i]} Data'],2))
TM = np.array(ListExtract(vars()[f'{files[i]} Data'],3))
n = np.array(ListExtract(vars()[f'{files[i]} Data'],6))
a = np.array(ListExtract(vars()[f'{files[i]} Data'],7)) / 1e7
kappa = np.array(ListExtract(vars()[f'{files[i]} Data'],10))


for i in range(len(files1)):
    files1[i] = files1[i].rstrip(".txt")
    
    vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    
    #This can be commented out with an alternet approach
    vars()[files1[i]][0] = vars()[files1[i]][0]
    
    vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
    vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],200,800)
    
    vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])

    s = F.Sub_nFinder(vars()['yP'+files1[i]])
    
x = xCalc(a,d)
phi = phiCalc(n,d,xP)
A = A1calc(n,s,kappa)
B = B1calc(n,s,kappa)
C = C1calc(n,s,kappa,phi)
D = D1calc(n,s,kappa)

T = ( (A*x) / (B - (C*x*np.cos(phi)) + (D * (x**2))) )

plt.figure(1)
plt.plot(xP,Traw, label = 'Actual')
plt.plot(xP,T, '--',label = 'Theoretical?')
plt.legend()