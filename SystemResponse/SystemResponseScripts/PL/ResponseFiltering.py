import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import copy
import random
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi
from whittaker_eilers import WhittakerSmoother
from sklearn.metrics import r2_score

def TrimSingle(Array,start,end):
    delPoints = np.where((Array < start) | (Array > end))[0]
    
    newArray = np.delete(Array,delPoints)
    
    return newArray


def Binning(array, n):
    return array[:(array.size // n) * n].reshape(-1, n).mean(axis=1)
    
def Man(Jonkles):
    ExpArr = np.floor(np.log10(abs(Jonkles)))
    
    Loc = np.where(ExpArr > (np.mean(ExpArr) + 3))[0]
    
    if Loc.size == 0:
        pass
    
    else:
        
        if Loc[0] == 0:
            Val = ExpArr[0] - ExpArr[1]
            Div = 10**Val
            Jonkles[0] = Jonkles[0] / Div
        
        if Loc[-1] == (len(ExpArr)-1):
            Val = ExpArr[-1] - ExpArr[-2]
            Div = 10**Val
            Jonkles[-1] = Jonkles[-1] / Div
    
        for i in range(len(Loc)):
            Val = ExpArr[Loc[i]] - np.floor(np.mean(ExpArr))
            Div = 10**Val
            Jonkles[Loc[i]] = Jonkles[Loc[i]] / Div
    
    
    return Jonkles

def RemoveZero(Array):
    
    if Array[0] == 0:
        Array[0] == np.where(Array != 0)[0][0]
    if Array[-1] == 0:
        Array[-1] == np.where(Array != 0)[0][-1]
    
    Locs = np.where(Array == 0)[0]
    for i in range(len(Locs)):
        Array[Locs[i]] = ((Array[(Locs[i]+1)] + Array[(Locs[i]-1)]) /2)
        
    return Array
        
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


def ConvertDependance(x,y):
    
    h = 6.62607015e-34
    c = 299792458
    x = x*1e-9
    
    xNew = ((h*c) / (x)) * 6.242e18
    yNew = y * (((x)**2)/(h*c))

    return xNew,yNew


def SpikeCorrect(xArray,yArray,xStart,xEnd,samples):
    Locs1 = np.where((xArray >= (xStart - samples)) & (xArray <= xStart))[0]
    Locs2 = np.where((xArray >= xEnd) & (xArray <= (xEnd + samples)))[0]
    Locs1 = np.append(Locs1,Locs2)  

    meanY = np.mean(yArray[Locs1])
    maxY = max(yArray[Locs1])
    minY = min(yArray[Locs1])
    
    Locs3 = np.where((xArray >= xStart) & (xArray <= xEnd))[0]
    TargetY = yArray[Locs3]
    
    for i in range(len(TargetY)):
        if TargetY[i] <= maxY and TargetY[i] >= minY:
            pass
        else:
            EndCon = 0
            while EndCon == 0:
                if TargetY[i] > maxY:
                    div = TargetY[i] / meanY
                    divNew = div + random.uniform(-0.25,0.25)
                    TargetY[i] = meanY * divNew
                elif TargetY[i] <  minY:
                    div = TargetY[i] / meanY
                    divNew = div + random.uniform(-0.25,0.25)
                    TargetY[i] = meanY * divNew                    
                else:
                    EndCon = 1
    
    newYarr = copy.deepcopy(yArray)
    
    for i in range(len(TargetY)):
        newYarr[Locs3[i]] = TargetY[i]
                
    
    return newYarr
    

def txtGrabber(scriptName,filePath,skipheader):
    
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))
    
    data = []
    
    #The next three loops remove the .txt from the files strings, this is for a cleaner title for the graph
    for i in range(len(fileNames)):
        fileNames[i] = fileNames[i][:-4]
        
        if skipheader == True:
            vars()[fileNames[i]] = np.loadtxt(open(path + "\\" + fileNames[i] + ".txt", "rb"), delimiter=",",skiprows = 2).T
        else:
            vars()[fileNames[i]] = np.loadtxt(open(path + "\\" + fileNames[i] + ".txt", "rb"), delimiter=",").T
            
        data.append(vars()[fileNames[i]])
    
    return fileNames, data
        

def fileGrabber(scriptName,filePath,suffix):
    
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))
    
    if len(fileNames) == 0:
        
        header = 'No Data'
        data = 0
        files = fileNames
    
    else:
    
        files = []
        
        for i in range(len(fileNames)):
            files.append(fileNames[i][:-len(suffix)])
            
        for i in range(len(files)):    
            rows = []
        
            with open(str(path)+"\\"+files[i]+ suffix,'r', encoding='utf-8-sig', newline='') as f:
                csvreader = csv.reader(f)
                header = next(csvreader)
            
                for row in csvreader:
                    rows.append(row)
            
            vars()[files[i]+'_Data_'+path[-1]] = rows
            
            for j in range(len(vars()[files[i]+'_Data_'+path[-1]])):
                for k in range(len(vars()[files[i]+'_Data_'+path[-1]][j])):
                        vars()[files[i]+'_Data_'+path[-1]][j][k] = float(vars()[files[i]+'_Data_'+path[-1]][j][k])
                    
        data = []
        for i in range(len(files)):
            data.append(vars()[files[i]+'_Data_'+path[-1]])
                    
    
    
    return path,dirs,files,data,header               



def fileGrabberCSVNH(scriptName,filePath):
    
    suffix = ".csv"
    
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))
    
    if len(fileNames) == 0:

        data = 0
        files = fileNames
    
    else:
    
        files = []
        
        for i in range(len(fileNames)):
            files.append(fileNames[i][:-len(suffix)])
            
        for i in range(len(files)):    
            rows = []
        
            with open(str(path)+"\\"+files[i]+ suffix,'r', encoding='utf-8-sig', newline='') as f:
                csvreader = csv.reader(f)
                
                for row in csvreader:
                    rows.append(row)
            
            vars()[files[i]+'_Data_'+path[-1]] = rows
            
            for j in range(len(vars()[files[i]+'_Data_'+path[-1]])):
                for k in range(len(vars()[files[i]+'_Data_'+path[-1]][j])):
                        vars()[files[i]+'_Data_'+path[-1]][j][k] = float(vars()[files[i]+'_Data_'+path[-1]][j][k])
                    
        data = []
        for i in range(len(files)):
            data.append(vars()[files[i]+'_Data_'+path[-1]])
                    
    
    
    return path,dirs,files,data 

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

def NoiseFilter(Val1,Val2,yArray):
    b, a = butter(Val1, Val2)
    yFiltered = filtfilt(b, a, yArray)
    
    return yFiltered

CLFileNames, CLData = txtGrabber('ResponseFiltering.py', '\\CLResponse', False)


pathQuanta1,dirsQuanta1,filesQuanta1,dataQuanta1  = fileGrabberCSVNH('ResponseFiltering.py', "\\LampInQuanta\\OceanOpticsLamp")
pathQuantaAvg,dirsQuantaAvg,filesQuantaAvg,dataQuantaAvg  = fileGrabberCSVNH('ResponseFiltering.py', "\\LampInQuanta\\OceanOpticsLamp\\Avg")
AndorFileNames, AndorData = txtGrabber('ResponseFiltering.py', '\\Andor', False)


for i in range(len(dataQuantaAvg)):
    vars()[f"DataSet {[i]} x"] = np.flip(np.array(ListExtract(dataQuantaAvg[i],0)))
    vars()[f"DataSet {[i]} y"] = np.flip(np.array(ListExtract(dataQuantaAvg[i],1)))
    

AvgDataSetY = np.zeros_like(vars()[f"DataSet {[0]} y"])

for i in range(len(vars()[f"DataSet {[0]} x"])):
    total = 0
    for j in range(len(dataQuantaAvg)):
        total += vars()[f"DataSet {[j]} y"][i]
        mean = (total / len(dataQuantaAvg))
        
    AvgDataSetY[i] = mean


LongTakeQuantaX = np.flip(np.array(ListExtract(dataQuanta1[0],0)))
LongTakeQuantaY = np.flip(np.array(ListExtract(dataQuanta1[0],1)))

AvgDataSetY = RemoveZero(AvgDataSetY)
LongTakeQuantaY = RemoveZero(LongTakeQuantaY)


NewYforXSysResp = np.interp(vars()[f"DataSet {[0]} x"],CLData[0][0],CLData[0][1])

LongTakeQuantaY = LongTakeQuantaY / max(LongTakeQuantaY)
AvgDataSetY = AvgDataSetY / max(AvgDataSetY)
NewYforXSysResp = NewYforXSysResp / max(NewYforXSysResp)

CorrectedData = AvgDataSetY / NewYforXSysResp

CorrectedData2 = LongTakeQuantaY / NewYforXSysResp

PLsetupX = AndorData[0][0]
PLsetupY = AndorData[0][1]

PLsetUPyNew = np.interp(vars()[f"DataSet {[0]} x"],AndorData[0][0],AndorData[0][1])

plt.plot(AndorData[0][0],AndorData[0][1])
plt.plot(vars()[f"DataSet {[0]} x"],PLsetUPyNew)

SystemResponse = PLsetUPyNew / CorrectedData
SystemResponse2 = PLsetUPyNew / CorrectedData2

SystemResponse = Man(SystemResponse)
SystemResponse2 = Man(SystemResponse2)

SystemResponseOriginal = copy.deepcopy(SystemResponse)
SystemResponseOriginal2 = copy.deepcopy(SystemResponse2)

plt.figure(6)
plt.plot(vars()[f"DataSet {[0]} x"],CorrectedData)
plt.plot(vars()[f"DataSet {[0]} x"],CorrectedData2)
plt.ylim(-1e-5,1e-5)
plt.ylim(-1,1)

TickSize = 16
LabelSize = 16
LegSize = 14
BinCount = 3
ValA = 3
ValB = 0.15


SystemResponse = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse,650,665,10)
SystemResponse2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse2,650,665,10)


SystemResponse = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse,650,665,10)
SystemResponse2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse2,650,665,10)

SystemResponse = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse,484,492,10)
SystemResponse2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponse2,484,492,10)

NewX = vars()[f"DataSet {[0]} x"][np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]
SystemResponse = SystemResponse[np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]
SystemResponse2 = SystemResponse2[np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]

lam = 3e5
smoothOrder = 2

ws = WhittakerSmoother(lmbda=lam, order=smoothOrder, data_length = len(SystemResponse))
yFilt11 = ws.smooth(SystemResponse)
yFilt12 = ws.smooth_optimal(SystemResponse, break_serial_correlation=False)
optimally_smoothed_series = yFilt12.get_optimal().get_smoothed()
optimal_lambda = yFilt12.get_optimal().get_lambda()

ws2 = WhittakerSmoother(lmbda=lam, order=smoothOrder, data_length = len(SystemResponse2))
yFilt13 = ws2.smooth(SystemResponse2)
yFilt14 = ws2.smooth_optimal(SystemResponse2, break_serial_correlation=False)
optimally_smoothed_series2 = yFilt14.get_optimal().get_smoothed()
optimal_lambda2 = yFilt14.get_optimal().get_lambda()



# SystemResponse = Binning(SystemResponse,BinCount)
# SystemResponse2 = Binning(SystemResponse2,BinCount)

# SystemResponse = NoiseFilter(ValA,ValB,SystemResponse)
# SystemResponse2 = NoiseFilter(ValA,ValB,SystemResponse2)

# xBinned = Binning(vars()[f"DataSet {[0]} x"],BinCount)

CheckFileName, CheckData = txtGrabber('ResponseFiltering.py', '\\CurrentlyUsedSystemResponse', False)

# NewX = vars()[f"DataSet {[0]} x"][np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]

# SystemResponse = SystemResponse[np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]
# SystemResponse2 = SystemResponse2[np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]

plt.figure(7)
plt.plot(vars()[f"DataSet {[0]} x"],SystemResponseOriginal)
plt.plot(NewX,SystemResponse)
plt.plot(NewX,yFilt11)
# plt.plot(vars()[f"DataSet {[0]} x"],SystemResponseOriginal2)
# plt.plot(CheckData[0][0],CheckData[0][1]/max(CheckData[0][1]))
plt.ylim(0,2e4)
plt.xlim(min(NewX),max(NewX))

# np.savetxt('ResponseOutputAvg.txt',np.c_[NewX,yFilt11], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLong.txt',np.c_[NewX,yFilt13], delimiter=",", fmt='%.3f')

# np.savetxt('ResponseOutputAvgUnfilt.txt',np.c_[NewX,SystemResponse], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLongUnfilt.txt',np.c_[NewX,SystemResponse2], delimiter=",", fmt='%.3f')


# SystemResponseOriginal = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponseOriginal,650,665,10)
# SystemResponseOriginal2 = SpikeCorrect(vars()[f"DataSet {[0]} x"],SystemResponseOriginal2,650,665,10)

# SaveNewX = vars()[f"DataSet {[0]} x"][np.where(vars()[f"DataSet {[0]} x"] >= 280)[0]]
# SaveNewY = np.interp(SaveNewX,xBinned,SystemResponse)
# SaveNewY2 = np.interp(SaveNewX,xBinned,SystemResponse2)

# SaveNewYUnfil = np.interp(SaveNewX,vars()[f"DataSet {[0]} x"],SystemResponseOriginal)
# SaveNewY2Unfil = np.interp(SaveNewX,vars()[f"DataSet {[0]} x"],SystemResponseOriginal2)

# np.savetxt('ResponseOutputAvg.txt',np.c_[SaveNewX,SaveNewY], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLong.txt',np.c_[SaveNewX,SaveNewY2], delimiter=",", fmt='%.3f')

# np.savetxt('ResponseOutputAvgUnfilt.txt',np.c_[SaveNewX,SaveNewYUnfil], delimiter=",", fmt='%.3f')
# np.savetxt('ResponseOutputLongUnfilt.txt',np.c_[SaveNewX,SaveNewY2Unfil], delimiter=",", fmt='%.3f')