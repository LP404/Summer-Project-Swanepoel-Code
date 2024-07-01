import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import nPlotterFunctions as nPF
import os
import SpectralResponseInterpArray as SRIA

def Normalise(x):
    
    xMax = max(x)
    xMin = min(x)
    xDiff = xMax - xMin
    
    xNorm = np.zeros_like(x)
    
    for i in range(len(x)):
        xNorm[i] = (x[i] - xMin) / xDiff
    
    
    
    return xNorm

def Trim(Array,start,end):
    delPoints = np.where((Array[0] < start) | (Array[0] > end))[0]
    
    
    newArray = np.array([np.delete(Array[0],delPoints),np.delete(Array[1],delPoints)])
    
    return newArray


def PlankLaw(T,x):
    h = 6.62607015e-34
    c = 299792458
    k = 1.380649e-23
    x = x*1e-9
    x = c/x
    
    A1 = (2*h*(x**3)) / (c**2)
    expo = ((h*x) / (k*T))
    A2 = 1 / (np.exp(expo) - 1)
    
    return A1*A2
    

GratingFileNames, GratingData = nPF.txtGrabber('SpectralResponse.py','\\SpectralResponse\\GratingComparison',False)
AndorFileNames, AndorData = nPF.txtGrabber('SpectralResponse.py', '\\SpectralResponse\\Andor', False)
OOFileNames, OOData = nPF.txtGrabber('SpectralResponse.py', '\\SpectralResponse\\OceanOptics', False)

TLFileNames, TLData = nPF.txtGrabber('SpectralResponse.py', '\\SpectralResponse\\TLamp', False)

CLFileNames, CLData = nPF.txtGrabber('SpectralResponse.py', '\\SpectralResponse\\SystemResponse\\CL', False)


pathA,dirsA,filesA,dataA  = nPF.fileGrabberCSVNH('SpectralResponse.py', '\\SpectralResponse\\SystemResponse\\CL\\SysteRes\\a')
pathB,dirsB,filesB,dataB  = nPF.fileGrabberCSVNH('SpectralResponse.py', '\\SpectralResponse\\SystemResponse\\CL\\SysteRes\\b')
pathBNew,dirsBNew,filesBNew,dataBNew  = nPF.fileGrabberCSVNH('SpectralResponse.py', '\\SpectralResponse\\SystemResponse\\CL\\SysteRes\\b\\New')
pathK,dirsK,filesK,dataK  = nPF.fileGrabberCSVNH('SpectralResponse.py', '\\SpectralResponse\\SystemResponse\\CL\\SysteRes\\k')


pathQuanta1,dirsQuanta1,filesQuanta1,dataQuanta1  = nPF.fileGrabberCSVNH('SpectralResponse.py', "\\SpectralResponse\\SystemResponse\\PL\\LampInQuanta\\OceanOpticsLamp")
pathQuantaAvg,dirsQuantaAvg,filesQuantaAvg,dataQuantaAvg  = nPF.fileGrabberCSVNH('SpectralResponse.py', "\\SpectralResponse\\SystemResponse\\PL\\LampInQuanta\\OceanOpticsLamp\\Avg")


# nPF.ListExtract(dataA[i],0)
# nPF.ListExtract(dataA[i],1)

h = 6.62607015e-34
c = 299792458

path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath('SpectralResponse.py')) + '\\PLTXT\\IWGO'))

for i in range(len(dirs)):
    vars()[f'files {i}'], vars()[f'data {i}'] = nPF.txtGrabber('SpectralResponse.py','\\PLTXT\\IWGO\\'+ dirs[i],True)

for i in range(len(dirs)):
    for j in range(len(vars()[f'files {i}'])):
        
        vars()[f'data {i}'][j][1] = vars()[f'data {i}'][j][1] * (((vars()[f'data {i}'][j][0]*1e-9)**2)/(h*c))
        vars()[f'data {i}'][j][0] = ((h*c) / (vars()[f'data {i}'][j][0]*1e-9)) * 6.242e18



xInterpMin = 290
xInterpMax = int(CLData[0][0][-1])
NoVals = 10001

xP = np.linspace(xInterpMin,xInterpMax,NoVals)



#77416 is the 400 @ 350 [0]
#77417 is the 400 @ 500 [1]
plt.figure(0)
plt.title('Values')
for i in range(len(GratingData)):
    plt.plot(GratingData[i][0]*1000,GratingData[i][1], label = GratingFileNames[i])
plt.legend()

plt.figure(1)
plt.title('Interp Values')
for i in range(len(GratingData)):
    vars()[f'yInterp {GratingFileNames[i]}'] = np.interp(xP,GratingData[i][0]*1000,GratingData[i][1])
    plt.plot(xP,vars()[f'yInterp {GratingFileNames[i]}'])

ScalarCorrection = vars()[f'yInterp {GratingFileNames[1]}'] / vars()[f'yInterp {GratingFileNames[0]}']
plt.figure(2)
plt.title('Interp Scalar Correction')
plt.plot(xP,1/ScalarCorrection)


plt.figure(3)
plt.title('77417 System Respone')
for i in range(len(CLData)):
    vars()[f'yInterp {CLFileNames[i]}'] = np.interp(xP,CLData[i][0],CLData[i][1])
    plt.plot(xP,vars()[f'yInterp {CLFileNames[i]}'])

Syst416 = vars()[f'yInterp {CLFileNames[i]}']/ScalarCorrection


plt.figure(4)
plt.title('Projected 77416 System Respone')
for i in range(len(CLData)):
    plt.plot(xP,Syst416)


CHIMParray = SRIA.ReturnArray()


clsArray = np.interp(CHIMParray,xP,Syst416)
cls2Array = np.interp(CHIMParray,CLData[0][0],CLData[0][1])

# np.savetxt('\\SpectralResponse\\SystemResponse\\PL\\ResponseTxt\\ResponseOutput.txt', np.c_[clsArray], fmt='%f')
# np.savetxt('\\SpectralResponse\\SystemResponse\\PL\\ResponseTxt\\ResponseOutput2.txt', np.c_[cls2Array], fmt='%f')

plt.figure(5)
plt.title('PL responses')
for i in range(len(AndorData)):
    plt.plot(AndorData[i][0],Normalise(AndorData[i][1]), label = 'Andor')
    plt.plot(OOData[i][0],Normalise(OOData[i][1]), label = 'Ocean Optics')
    plt.legend()

#xStart2 = AndorData[0][0][0]
xStart2 = 244
xEnd2 = AndorData[0][0][-1]
NoXvals2 = 10001

xP2 = np.linspace(xStart2,xEnd2,NoXvals2)

AndorInterp = np.interp(xP2,AndorData[0][0],Normalise(AndorData[0][1]))
OOInterp = np.interp(xP2,OOData[0][0],Normalise(OOData[0][1]))

plt.figure(6)
plt.title('PL Interp')
plt.plot(xP2,AndorInterp, label = 'Andor')
plt.plot(xP2,OOInterp, label = 'Ocean Optics')
plt.legend()

Response = AndorInterp/OOInterp

plt.figure(7)
plt.title('PL Delta')
plt.plot(xP2,Response)
plt.legend()

plt.figure(8)
plt.xlim(200,700)
plt.title('Tungsten Lamp')
plt.plot(TLData[0][0],Normalise(TLData[0][1]), label = 'Experemental')
plt.plot(TLData[0][0],Normalise(PlankLaw(2400,TLData[0][0])), label = 'Theoretical')


TLInterp = np.interp(xP2,TLData[0][0],Normalise(TLData[0][1]))
BBInterp = np.interp(xP2,TLData[0][0],Normalise(PlankLaw(2400,TLData[0][0])))

Response2 = BBInterp/TLInterp

plt.figure(9)
plt.xlim(300,AndorData[0][0][-1])
plt.ylim(0,1)
plt.title('Interp')
plt.plot(xP2,TLInterp, label = 'Experemental')
plt.plot(xP2,BBInterp, label = 'Theoretical')
plt.plot(xP2,Response2, label = 'Delta')
plt.legend()

plt.figure(10)
plt.title('Response Convolution')
plt.plot(xP2,(Response2*Response), label = 'Delta')
plt.legend()