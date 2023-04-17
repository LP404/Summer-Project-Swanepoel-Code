import numpy as np

def PerToAbs(PerUncert,Val):
    AbsUncert = (PerUncert/100)*Val
    return AbsUncert

def AbsToPer(AbsUncert,Val):
    PerUncert = (AbsUncert/Val)*100
    return PerUncert

#Used to add uncetanties together for either absolutes for add/subtract or percentages for multuplt/deivide
def LinearCombine(Uncert1,Uncert2):
    CombiUncert = abs(Uncert1) + abs(Uncert2)
    return CombiUncert

def MultiDivUncert(Uncert1,Uncert2,Val1,Val2):
    CombiUncert = ((abs(Uncert1)/Val1)*100) + ((abs(Uncert2)/Val2)*100)
    return CombiUncert

def ExpoUncert(Uncert,Val,Expo):
    NewUncert = Expo*((abs(Uncert)/Val)*100)
    return NewUncert

def ScalingUncert(Uncert,Scalar):
    return (Uncert*Scalar)

TUncert = 0.5/100
TM = 0.80
Tm = 0.75
n = 1.97


Lam1 = 450
Lam2 = 475

n1 = 1.98
n2 = 1.99

LamUncert = 0.3

n1Uncert = 0.017
n2Uncert = 0.0191

d = 1600

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