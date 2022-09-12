import matplotlib.pyplot as plt
import csv
import os
import numpy as np
from scipy.optimize import curve_fit , leastsq
from scipy.stats import linregress


def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Plotter.py')) + '\\CSVout'))
path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('Plotter.py')) + '\\CSVout\\Additional_Info'))

suffix = ".csv"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  

    rows = []
    
    with open(str(path)+"\\"+files[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header = next(csvreader)
        
        for row in csvreader:
            rows.append(row)
    
    vars()[files[i]+'Data'] = rows

    for j in range(len(vars()[files[i]+'Data'])):
        
        for k in range(len(vars()[files[i]+'Data'][j])):
        
            if vars()[files[i]+'Data'][j][k] == "yes":
                vars()[files[i]+'Data'][j][k] = True
                
            elif vars()[files[i]+'Data'][j][k] == "no":
                vars()[files[i]+'Data'][j][k] = False
                
            else:
                vars()[files[i]+'Data'][j][k] = float(vars()[files[i]+'Data'][j][k])

    vars()[files[i]+'Delte_Rows'] = []
    
    for j in range(len(vars()[files[i]+'Data'])):

        if vars()[files[i]+'Data'][j][10] == True:    
            vars()[files[i]+'Delte_Rows'].append(j)
     
    vars()[files[i]+'Delte_Rows'] = sorted(vars()[files[i]+'Delte_Rows'],reverse=True)
    
    for j in range(len(vars()[files[i]+'Delte_Rows'])):
        del vars()[files[i]+'Data'][vars()[files[i]+'Delte_Rows'][j]]



for i in range(len(files1)):
    files1[i] = files1[i][:-len(suffix)]  

    rows1 = []
    
    with open(str(path1)+"\\"+files1[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f1:
        csvreader1 = csv.reader(f1)
        header1 = next(csvreader1)
    
        for row1 in csvreader1:
            rows1.append(row1)
    
    vars()[files1[i]+'_Data'] = rows1
    
    for j in range(len(vars()[files1[i]+'_Data'])):
        for k in range(len(vars()[files1[i]+'_Data'][j])):
                vars()[files1[i]+'_Data'][j][k] = float(vars()[files1[i]+'_Data'][j][k])
                
for i in range(len(files)):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.figure(i)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='darkgray', linestyle='--')
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],3), label = files[i]+" n1")
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],7), label = files[i]+" n2")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Refractive Index")
    plt.legend()
    
  
    n = np.array(ListExtract(vars()[files[i]+'Data'],7))
    Lam = np.array(ListExtract(vars()[files[i]+'Data'],0))
    m1 = ListExtract(vars()[files[i]+'Data'],5)[len(ListExtract(vars()[files[i]+'Data'],5))-1]
    y =  np.arange(0,len(n),1) / 2
    x = n/Lam
    vars()[files[i]+'FitVals'] = linregress(x,y)
    fx = np.linspace(0,2*max(x),10001)
    fy = vars()[files[i]+'FitVals'][0] * fx + vars()[files[i]+'FitVals'][1]
    plt.figure(100+i)
    plt.title(files[i])
    plt.ylabel("l/2")
    plt.xlabel('n/λ')
    plt.scatter(x,y)
    plt.plot(fx,fy)
    

    plt.figure(200+i)
    plt.title(files[i]+ 'Extinction Coeffiecent')
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],8))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("κ")
    
 
    plt.figure(300+i)
    plt.title(files[i]+ 'Absorbtion Coeffiecent')
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],9))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("α (cm^-1)")