import pandas as pd
import numpy as np
import os
import csv
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib as mpl
import nPlotterFunctions as nPF
from scipy.interpolate import CubicSpline, PchipInterpolator, Akima1DInterpolator

CLFileNames, CLData = nPF.txtGrabber('CLcorrector.py', '\\SystemResponse\\SystemResponseScripts\\CL', False)
xResp,yResp = nPF.ConvertDependance(CLData[0][0],CLData[0][1])

xResp = np.flip(xResp)
yResp = np.flip(yResp)

scriptName = 'CLcorrector.py'
filePath = '\\CLdata\\Input\\b'

# path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))

# df1 = pd.read_excel(path + f"\{fileNames[0]}", header=None).T.to_numpy()

def xlsxGrabber(scriptName,filePath,skipheader):
    
    path, dirs, fileNames = next(os.walk(os.path.dirname(os.path.realpath(scriptName)) + filePath))
    
    data = []
    

    for i in range(len(fileNames)):
        #Removes the .xlsx extentio
        exten = fileNames[i][fileNames[i].find('.'):]
        
        fileNames[i] = fileNames[i][0:fileNames[i].find('.')]
        
        if skipheader == True:
            vars()[fileNames[i]] = pd.read_excel(path + f"\{fileNames[i]}" + exten, header=None).T.to_numpy()
        else:
            vars()[fileNames[i]] = pd.read_excel(path + f"\{fileNames[i]}" + exten).T.to_numpy()
            
        data.append(vars()[fileNames[i]])
    
    return fileNames, data


files,data = xlsxGrabber(scriptName,filePath,True)

filesSorted = natsorted(files)
NewData = [x for _, x in natsorted(zip(files, data))]

GraphLabels = []
for i in range(len(filesSorted)):
    GraphLabels.append(filesSorted[i][:filesSorted[i].find('_')])
    
    
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['xkcd:blue', 'xkcd:forest green', 'xkcd:dark red', 'xkcd:cyan', 'xkcd:magenta', 'xkcd:black', 'xkcd:purple', 'xkcd:brown', 'xkcd:orange', 'xkcd:teal', 'xkcd:coral', 'xkcd:light blue', 'xkcd:lime', 'xkcd:turquoise', 'xkcd:dark green', 'xkcd:tan', 'xkcd:salmon', 'xkcd:gold', 'xkcd:grey', 'xkcd:pink', 'xkcd:apricot','xkcd:olive','xkcd:navy','xkcd:lavender','xkcd:beige','xkcd:teal'])




plt.figure(0,figsize = (16,10))
for i in range(len(filesSorted)):
    if np.random.default_rng().random() > 0.5:
       varia = 'dashed'
    else:
       varia = 'solid'    
    plt.semilogy(NewData[i][0],NewData[i][1], label = GraphLabels[i], linestyle = varia)
plt.legend()
plt.xlim(min(NewData[0][0]),4.5)



yResponse = Akima1DInterpolator(xResp,yResp)(NewData[0][0])

plt.figure(1,figsize = (16,10))
for i in range(len(filesSorted)):
    if np.random.default_rng().random() > 0.5:
       varia = 'dashed'
    else:
       varia = 'solid'    
    plt.semilogy(NewData[i][0],NewData[i][1]/yResponse, label = GraphLabels[i], linestyle = varia)
plt.legend()
plt.xlim(min(NewData[0][0]),4.5)

#a/b
# plt.ylim(0,1e-16)
#k
