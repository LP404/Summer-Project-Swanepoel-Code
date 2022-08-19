import matplotlib.pyplot as plt
import csv
import os
from operator import itemgetter

# def ListSort(List):
#     return(sorted(List, map(itemgetter(0), List ))) 

#Swanepoel Range Control
SRC = 2

def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

path, dirs, files = next(os.walk(os.path.dirname(os.path.realpath('Code.py')) + '\\CSVout'))

suffix = ".csv"
for i in range(len(files)):
    files[i] = files[i][:-len(suffix)]  
    
 
for i in range(len(files)):  
    rows = []
    with open(str(path)+"\\"+files[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header = next(csvreader)
        for row in csvreader:
            rows.append(row)
    
    vars()[files[i]+'Data'] = rows
    
for i in range(len(files)):  
    for j in range(len(vars()[files[i]+'Data'])):
        for k in range(len(vars()[files[i]+'Data'][j])):
            if vars()[files[i]+'Data'][j][k] == "-":
                pass
            else:
                vars()[files[i]+'Data'][j][k] = float(vars()[files[i]+'Data'][j][k])
                
for i in range(len(files)): 
    vars()[files[i]+'Data'] = ListSort(vars()[files[i]+'Data'])
    del vars()[files[i]+'Data'][0:SRC]
    del vars()[files[i]+'Data'][len(vars()[files[i]+'Data'])-SRC:len(vars()[files[i]+'Data'])]
    
    
for i in range(len(files)):
    # plt.figure(i,figsize=(29.7/2.54,21.0/2.54), dpi=600)
    plt.figure(i)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='darkgray', linestyle='--')
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],3), label = files[i]+"_RI1")
    plt.scatter(ListExtract(vars()[files[i]+'Data'],0),ListExtract(vars()[files[i]+'Data'],7), label = files[i]+"_RI2")
    plt.title(files[i])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Refractive Index")
    plt.legend()