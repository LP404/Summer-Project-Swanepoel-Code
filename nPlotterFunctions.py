import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from scipy import stats
from scipy.optimize import curve_fit
import math
import Functions as F
from sklearn.metrics import r2_score

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]

# #GrabsAllCSVs in a Folder
def fileGrabberTrans(scriptName,filePath,suffix):
    
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
                    try:
                        vars()[files[i]+'_Data_'+path[-1]][j][k] = float(vars()[files[i]+'_Data_'+path[-1]][j][k])
                    except ValueError:
                        pass
        
        
        for i in range(len(files)): 
            MaxLim = len(vars()[files[i]+'_Data_'+path[-1]])
            j = 0
            while j < MaxLim:
                if vars()[files[i]+'_Data_'+path[-1]][j][18] == 'yes':
                    vars()[files[i]+'_Data_'+path[-1]].pop(j)
                    MaxLim = len(vars()[files[i]+'_Data_'+path[-1]])
                else:
                    j +=1
            
        data = []
        for i in range(len(files)):
            data.append(vars()[files[i]+'_Data_'+path[-1]])
                    
    
    
    return path,dirs,files,data,header
                    
        
def datGrabber(datName,filepath):
    data = np.loadtxt(filepath + datName + '.dat')
    data = data.T
    
    if len(data) == 10:
        header = ['x','y','rho','active','peak','peak convolution','x-correction','residuals','absolute residuals','weigthed residuals']
    else:
        diff = len(data) - 10
        header = ['x','y','rho','active','peak 1','peak convolution','x-correction','residuals','absolute residuals','weigthed residuals',diff + 1]
        for i in range(diff):
            header.insert(i+5,f'peak {i+2}')
            
    return header,data

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


def CauchyFinder(files,data,xP,fileIndPos,fileIndPos1,fileIndPos2): 
    

    vars()['CauchFindY1 '+files] = np.array(ListExtract(data,fileIndPos))
    vars()['CauchFindY2 '+files] = np.array(ListExtract(data,fileIndPos1))
    vars()['CauchFindX '+files] = 1/((np.array(ListExtract(data,fileIndPos2)))**2)
    
    n1Cauch = [vars()['CauchFindY1 '+files],vars()['CauchFindX '+files]]
    n2Cauch = [vars()['CauchFindY2 '+files],vars()['CauchFindX '+files]]
    
    vars()[f"mn1-{files}"], vars()[f"interceptn1-{files}"], vars()[f"r value n1-{files}"], vars()["p value n1-{files}"], vars()[f"std errn1-{files}"] = stats.linregress(vars()['CauchFindX '+files], vars()['CauchFindY1 '+files])
    vars()[f"mn2-{files}"], vars()[f"interceptn2-{files}"], vars()[f"r value n2-{files}"], vars()["p value n2-{files}"], vars()[f"std errn2-{files}"] = stats.linregress(vars()['CauchFindX '+files], vars()['CauchFindY2 '+files])
       
    n1CoefList = [vars()[f"mn1-{files}"], vars()[f"interceptn1-{files}"], vars()[f"r value n1-{files}"], vars()["p value n1-{files}"], vars()[f"std errn1-{files}"]]
    n2CoefList = [vars()[f"mn2-{files}"], vars()[f"interceptn2-{files}"], vars()[f"r value n2-{files}"], vars()["p value n2-{files}"], vars()[f"std errn2-{files}"]]
       
    vars()[f"yFitn1 {files}"] = vars()[f"mn1-{files}"] / (xP**2) + vars()[f"interceptn1-{files}"]
    vars()[f"yFitn2 {files}"] = vars()[f"mn2-{files}"] / (xP**2) + vars()[f"interceptn2-{files}"]
    
    n1yFit = vars()[f"yFitn1 {files}"]
    n2yFit = vars()[f"yFitn2 {files}"]
    
    return n1yFit, n2yFit, n1Cauch, n2Cauch, n1CoefList, n2CoefList
         
         
         