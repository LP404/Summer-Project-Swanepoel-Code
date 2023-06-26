import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import os
import csv


def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

path1, dirs1, files1 = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\cauchFit'))

suffix = ".csv"
for i in range(len(files1)):
    files1[i] = files1[i][:-len(suffix)]  
    
    rows1 = []

    with open(str(path1)+"\\"+files1[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header1 = next(csvreader)
    
        for row1 in csvreader:
            rows1.append(row1)
    
    vars()[files1[i]+'_Data'] = rows1
    
    for j in range(len(vars()[files1[i]+'_Data'])):
        for k in range(len(vars()[files1[i]+'_Data'][j])):
                vars()[files1[i]+'_Data'][j][k] = float(vars()[files1[i]+'_Data'][j][k])
                
                
x = np.array(ListExtract(vars()[files1[i]+'_Data'],0))   
y =  np.array(ListExtract(vars()[files1[i]+'_Data'],2))       
xP = np.arange(190,801,1)
xVal = np.array([])
yVal = np.array([])

for i in range(len(files1)):
    for k in range(len(xP)):
    
        Loc = np.where(x == F.FindNearestVal(x,xP[k]))[0][0]
        xVal = np.append(xVal,np.around(x[Loc],0))
        yVal = np.append(yVal,y[Loc])
        
        
    header = ['λ (nm)','n2 fit']

    with open('CSVout/'+files1[i]+ '_FitSimplified.csv','w', encoding='utf-8-sig', newline='') as f2:
         
        writer = csv.writer(f2)
        writer.writerow(header)
        
        for l in range(len(xVal)):
            data = [xVal[l],yVal[l]]
            writer.writerow(data)