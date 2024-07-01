import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import os
import csv
import nPlotterFunctions as nPF

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

def xCalc (a,d):
    return np.exp((-1*a*d))

def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]


def ktoa(k,lam):
    return (4*np.pi * k)/lam


path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\SubstrateData'))
pathN,dirsN,filesN,dataN,headerN  = nPF.fileGrabber('TheoreticalTransmittanceTester.py','\\n','.csv')
pathK,dirsK,filesK,dataK,headerK  = nPF.fileGrabber('TheoreticalTransmittanceTester.py','\\k','.csv')


for i in range(len(filesN)):
    vars()[f'{filesN[i]} xN'] = np.array(ListExtract(dataN[i],0))
    vars()[f'{filesN[i]} yN'] = np.array(ListExtract(dataN[i],1))

for i in range(len(filesN)):
    if max(vars()[f'{filesN[i]} xN']) < 200:
        vars()[f'{filesN[i]} xN'] = vars()[f'{filesN[i]} xN'] * 1000

for i in range(len(filesK)):
    vars()[f'{filesK[i]} xK'] = np.array(ListExtract(dataK[i],0))
    vars()[f'{filesK[i]} yK'] = np.array(ListExtract(dataK[i],1))

for i in range(len(filesK)):
    if max(vars()[f'{filesK[i]} xK']) < 200:
        vars()[f'{filesK[i]} xK'] = vars()[f'{filesK[i]} xK'] * 1000


for i in range(len(files1)):
    files1[i] = files1[i].rstrip(".txt")
    
    vars()[files1[i]] = np.loadtxt(open(path1 + "\\" + files1[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
    
    #This can be commented out with an alternet approach
    vars()[files1[i]][0] = vars()[files1[i]][0]
    
    vars()[files1[i]][1] = vars()[files1[i]][1] / 100
    
    vars()[files1[i]+'T'] = F.Trim(vars()[files1[i]],200,800)
    
    xP = np.linspace(200,1000,100001)
    
    vars()['yP'+files1[i]] = np.interp(xP,vars()[files1[i]+'T'][0],vars()[files1[i]+'T'][1])

    s = F.Sub_nFinder(vars()['yP'+files1[i]])

#Multipling 1e7 converts from nm^-1 to cm^-1
for i in range(len(filesK)):
    vars()[f'{filesK[i]} yA'] = ktoa(vars()[f'{filesK[i]} yK'],vars()[f'{filesK[i]} xK']) * 1e7

for i in range(len(filesK)):
    vars()[f'newA {i}'] = np.interp(xP,vars()[f'{filesK[i]} xK'], vars()[f'{filesK[i]} yA'])
    vars()[f'newK {i}'] = np.interp(xP,vars()[f'{filesK[i]} xK'], vars()[f'{filesK[i]} yK'])
    vars()[f'newN {i}'] = np.interp(xP,vars()[f'{filesN[i]} xN'], vars()[f'{filesN[i]} yN'])

d = np.arange(100,2100,100)

for j in range(len(d)):
    for i in range(len(filesK)):
        vars()[f'x {i}_{d[j]}'] = xCalc(vars()[f'newA {i}'],d[j])
        vars()[f'phi {i}_{d[j]}'] = phiCalc(vars()[f'newN {i}'],d[j],xP)
        vars()[f'A {i}_{d[j]}'] = A1calc(vars()[f'newN {i}'],s,vars()[f'newK {i}'])
        vars()[f'B {i}_{d[j]}'] = B1calc(vars()[f'newN {i}'],s,vars()[f'newK {i}'])
        vars()[f'C {i}_{d[j]}'] = C1calc(vars()[f'newN {i}'],s,vars()[f'newK {i}'],vars()[f'phi {i}_{d[j]}'])
        vars()[f'D {i}_{d[j]}'] = D1calc(vars()[f'newN {i}'],s,vars()[f'newK {i}'])

for j in range(len(d)):
    for i in range(len(filesK)):
            vars()[f'T {i}_{d[j]}'] = ( (vars()[f'A {i}_{d[j]}']*vars()[f'x {i}_{d[j]}']) / (vars()[f'B {i}_{d[j]}'] - (vars()[f'C {i}_{d[j]}']*vars()[f'x {i}_{d[j]}']*np.cos(vars()[f'phi {i}_{d[j]}'])) + (vars()[f'D {i}_{d[j]}'] * (vars()[f'x {i}_{d[j]}']**2))) )

for j in range(len(d)):
    for i in range(len(filesK)):
        plt.figure(i + 20**j)
        plt.plot(xP, vars()[f'T {i}_{d[j]}'],label = f'Theoretical {filesK[i]} and {d[j]}nm ')
        plt.legend()
        
for j in range(len(d)):
    for i in range(len(filesK)):
            vars()[f'T {i}_{d[j]}'] = ( (vars()[f'A {i}_{d[j]}']*vars()[f'x {i}_{d[j]}']) / (vars()[f'B {i}_{d[j]}'] - (vars()[f'C {i}_{d[j]}']*vars()[f'x {i}_{d[j]}']*np.cos(vars()[f'phi {i}_{d[j]}'])) + (vars()[f'D {i}_{d[j]}'] * (vars()[f'x {i}_{d[j]}']**2))) )
            np.savetxt(f'T {i}_{d[j]}.txt',np.column_stack((xP,vars()[f'T {i}_{d[j]}']*100)),delimiter=',')