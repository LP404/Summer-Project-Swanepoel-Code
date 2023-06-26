import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from scipy import stats
from scipy.optimize import curve_fit
import math
import Functions as F
from sklearn.metrics import r2_score
import nPlotterFunctions as nPF

def CTI(value):
    return value/2.54

def ListSort(List):
    return(sorted(List, key = lambda x: x[0]))   

def ListExtract(List,IndVal):
    return [Index[IndVal] for Index in List]


def choose_subplot_dimensions(k):
    if k < 4:
        return k, 1
    elif k < 11:
        return math.ceil(k/2), 2
    else:
        # I've chosen to have a maximum of 3 columns
        return math.ceil(k/3), 3


def sellmerier(lam,b1,b2,b3,c1,c2,c3):
    return np.sqrt(( 1 + ((b1*lam**2)/(lam**2 - c1)) + ((b2*lam**2)/(lam**2 - c2)) + ((b3*lam**2)/(lam**2 - c3))))

    
plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))




DeltaLam = 0.3

path, dirs, files, data, header = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\α','.csv')
path1, dirs1, files1, data1, header1 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\β','.csv')
path2, dirs2, files2, data2, header2 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\κ','.csv')

path3, dirs3, files3 = next(os.walk(os.path.dirname(os.path.realpath('BandGapCalc.py')) + '\\CSVout\\Absorption'))

suffix = ".csv"

for i in range(len(files3)):
    files3[i] = files3[i][:-len(suffix)]  
    
    rows1 = []

    with open(str(path3)+"\\"+files3[i]+ '.csv','r', encoding='utf-8-sig', newline='') as f:
        csvreader = csv.reader(f)
        header1 = next(csvreader)
    
        for row1 in csvreader:
            rows1.append(row1)
    
    vars()[files3[i]+'_Data'] = rows1
    
    for j in range(len(vars()[files3[i]+'_Data'])):
        for k in range(len(vars()[files3[i]+'_Data'][j])):
                vars()[files3[i]+'_Data'][j][k] = float(vars()[files3[i]+'_Data'][j][k])


#Ind 0,5,6 n1
#Ind 0,13,14


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','red','black','gold','orange','salmon','silver','pink','teal','purple']

xP = np.linspace(190,800, 10001)

if data != 0:
    for i in range(len(files)):
        vars()[f"n1yFit {files[i]}"], vars()[f"n2yFit {files[i]}"], vars()[f"n1Cauch {files[i]}"], vars()[f"n2Cauch {files[i]}"], vars()[f"n1CoefList {files[i]}"], vars()[f"n2CoefList {files[i]}"] = nPF.CauchyFinder(files[i],data[i],xP,5,13,0)
    #     plt.figure(i + 5*len(files))
    #     plt.minorticks_on()
    #     plt.grid(which='major', color='k', linestyle='-')
    #     plt.grid(which='minor', color='darkgray', linestyle='--')
    #     plt.plot(xP, vars()[f"n1yFit {files[i]}"], '--',color = LineColour[i],label = f'n1 = {np.around(vars()[f"n1CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n1CoefList {files[i]}"][0],3):.2e}/λ^2')
    #     plt.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
    #     plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],5),ListExtract(data[i],6),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
    #     plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
    #     plt.title(LabList[0] + 'Refrative Index')
    #     plt.ylabel('n')
    #     plt.xlabel('λ (nm)')
    #     plt.legend()   
    #     plt.plot()
        
    num = len(files)
    plt.subplots_adjust(hspace=1,wspace=1)
    SubRow,SubCol = choose_subplot_dimensions(num)
    plt.subplots(figsize=(15,12),sharex=True)
    plt.suptitle(LabList[0] + 'Refrative Index')
    for i in range(len(files)):

        # add a new subplot iteratively
        plt.minorticks_on()
        plt.grid(which='major', color='k', linestyle='-')
        plt.grid(which='minor', color='darkgray', linestyle='--')
        ax = plt.subplot(SubRow, SubCol, i + 1)

        # filter df and plot ticker on the new subplot axis
        ax.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        ax.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
        ax.legend()
        ax.title.set_text(f'n2 for {files[i]}')

    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.savefig('AlphaN2.png', dpi = 600)
    
    
    for i in range(len(files)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files[i]}"][k],vars()[f"n2yFit {files[i]}"][k]]
                
                writer.writerow(dataWrite)
     
    
else:
    pass

if data1 != 0:
    for i in range(len(files1)):
        vars()[f"n1yFit {files1[i]}"], vars()[f"n2yFit {files1[i]}"], vars()[f"n1Cauch {files1[i]}"], vars()[f"n2Cauch {files1[i]}"], vars()[f"n1CoefList {files1[i]}"], vars()[f"n2CoefList {files1[i]}"] = nPF.CauchyFinder(files1[i],data1[i],xP,5,13,0)
        # plt.figure(i + 14*len(files))
        # plt.minorticks_on()
        # plt.grid(which='major', color='k', linestyle='-')
        # plt.grid(which='minor', color='darkgray', linestyle='--')
        # plt.plot(xP, vars()[f"n1yFit {files1[i]}"], '--',color = LineColour[i],label = f'n1 = {np.around(vars()[f"n1CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n1CoefList {files1[i]}"][0],3):.2e}/λ^2')
        # plt.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
        # plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],5),ListExtract(data2[i],6),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        # plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        # plt.title(LabList[1] + 'Refrative Index')
        # plt.ylabel('n')
        # plt.xlabel('λ (nm)')
        # plt.legend()   
        # plt.plot()
        
    num1 = len(files1)
    plt.figure(2)
    plt.subplots_adjust(hspace=1,wspace=1)
    SubRow,SubCol = choose_subplot_dimensions(num1)
    plt.subplots(figsize=(15,12),sharex=True)
    plt.suptitle(LabList[1] + 'Refrative Index')
    for i in range(len(files1)):

        # add a new subplot iteratively
        plt.minorticks_on()
        plt.grid(which='major', color='k', linestyle='-')
        plt.grid(which='minor', color='darkgray', linestyle='--')
        ax = plt.subplot(SubRow, SubCol, i + 1)

        # filter df and plot ticker on the new subplot axis
        ax.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        ax.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
        ax.legend()
        ax.title.set_text(f'n2 for {files1[i]}')

    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.savefig('BetaN2.png', dpi = 600)
    
    
    for i in range(len(files1)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files1[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files1[i]}"][k],vars()[f"n2yFit {files1[i]}"][k]]
                
                writer.writerow(dataWrite)

else:
    pass

if data2 != 0:
    for i in range(len(files2)):
        vars()[f"n1yFit {files2[i]}"], vars()[f"n2yFit {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2CoefList {files2[i]}"], vars()[f"n2CoefList {files2[i]}"] = nPF.CauchyFinder(files2[i],data2[i],xP,5,13,0)
        # plt.figure(i + 17*len(files))
        # plt.minorticks_on()
        # plt.grid(which='major', color='k', linestyle='-')
        # plt.grid(which='minor', color='darkgray', linestyle='--')
        # plt.plot(xP, vars()[f"n1yFit {files2[i]}"], '--',color = LineColour[i],label = f'n1 = {np.around(vars()[f"n1CoefList {files2[i]}"][1],3)} + {np.around(vars()[f"n1CoefList {files2[i]}"][0],3):.2e}/λ^2')
        # plt.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files2[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files2[i]}"][0],3):.2e}/λ^2')
        # plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],5),ListExtract(data2[i],6),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        # plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        # plt.title(LabList[1] + 'Refrative Index')
        # plt.ylabel('n')
        # plt.xlabel('λ (nm)')
        # plt.legend()   
        # plt.plot()
        
    # num2 = len(files2)
    plt.figure(3)
    plt.subplots_adjust(hspace=1,wspace=1)
    SubRow,SubCol = choose_subplot_dimensions(num1)
    plt.subplots(figsize=(15,12),sharex=True)
    plt.suptitle(LabList[2] + 'Refrative Index')
    for i in range(len(files2)):

        # add a new subplot iteratively
        plt.minorticks_on()
        plt.grid(which='major', color='k', linestyle='-')
        plt.grid(which='minor', color='darkgray', linestyle='--')
        ax = plt.subplot(SubRow, SubCol, i + 1)

        # filter df and plot ticker on the new subplot axis
        ax.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        ax.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+1],label = f'n2 = {np.around(vars()[f"n2CoefList {files2[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files2[i]}"][0],3):.2e}/λ^2')
        ax.legend()
        ax.title.set_text(f'n2 for {files2[i]}')

    plt.minorticks_on()
    plt.grid(which='major', color='k', linestyle='-')
    plt.grid(which='minor', color='darkgray', linestyle='--')
    plt.savefig('BetaN2.png', dpi = 600)
    
    for i in range(len(files2)):
        header = ['λ (nm)','n1 fit','n2 fit']
        with open('CSVout/cauchFit/'+files2[i]+ '_CauchFit.csv','w', encoding='utf-8-sig', newline='') as f:
            
            writer = csv.writer(f)
            
            writer.writerow(header)
        
            for k in range(len(xP)):
                dataWrite = [xP[k],vars()[f"n1yFit {files2[i]}"][k],vars()[f"n2yFit {files2[i]}"][k]]
                
                writer.writerow(dataWrite)
else:
    pass

LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']

plt.figure(40)
plt.minorticks_on()
plt.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
plt.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
plt.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
plt.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
plt.title('n values for different Ga₂O₃ polymorphs')
plt.ylabel('n (arb. units)')
plt.xlabel(r'$\lambda$ (nm)')
for i in range(len(files)):
        plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        plt.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
for i in range(len(files1)):
        plt.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        plt.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
plt.legend()

    
# for i in range(len(files)):
#     Loc = np.where((xP >= 380) & (xP<= 700))[0]
    
#     vars()[f'SfitY {files[i]}'] = vars()[f"yFitn2 {files[i]}"][Loc]
#     vars()[f'SfitX {files[i]}'] = xP[Loc]
    
#     Loc2 = np.where((np.array(ListExtract(vars()[files[i]+'_Data'],0)) < 380) | (np.array(ListExtract(vars()[files[i]+'_Data'],0)) > 700))[0]
    
#     vars()[f'SfitY {files[i]}'] = np.append(vars()[f'SfitY {files[i]}'],np.array(ListExtract(vars()[files[i]+'_Data'],13))[Loc2][0])
#     vars()[f'SfitX {files[i]}'] = np.append(vars()[f'SfitX {files[i]}'],np.array(ListExtract(vars()[files[i]+'_Data'],0))[Loc2][0])   

#     LambdaInds = vars()[f'SfitX {files[i]}'].argsort()
    
#     vars()[f'SfitY {files[i]}'] = vars()[f'SfitY {files[i]}'][LambdaInds[::1]]
#     vars()[f'SfitX {files[i]}'] = vars()[f'SfitX {files[i]}'][LambdaInds[::1]]
    
#     #Guess based off of sapphire values
#     guess = [1.99,0.15,2.84,13280,32200,6.25e9]
    
#     c,cov = curve_fit(sellmerier,vars()[f'SfitX {files[i]}'],vars()[f'SfitY {files[i]}'],guess,maxfev=5000000)

    
#     GuessY = sellmerier(vars()[f'SfitX {files[i]}'],guess[0],guess[1],guess[2],guess[3],guess[4],guess[5])
#     CalcY = sellmerier(vars()[f'SfitX {files[i]}'],c[0],c[1],c[2],c[3],c[4],c[5])
#     FitY = sellmerier(xP,c[0],c[1],c[2],c[3],c[4],c[5])
    

# xP2 = np.linspace(190,1500,100001)
# yFit2extrap = vars()[f"mn2-{files[i]}"] / (xP2**2) + vars()[f"interceptn2-{files[i]}"]
# FitY2 = sellmerier(xP2,c[0],c[1],c[2],c[3],c[4],c[5])

# yVal = np.array(ListExtract(vars()[files[i]+'_Data'],13))
# xVal = np.array(ListExtract(vars()[files[i]+'_Data'],0))

# yCheckC = np.array([])
# yCheckS = np.array([])

# for i in range(len(xVal)):
#     Loc = np.where(xP2 == F.FindNearestVal(xP2,xVal[i]))[0][0]
#     yCheckC = np.append(yCheckC,yFit2extrap[Loc])
#     yCheckS = np.append(yCheckS,FitY2[Loc])

# rScoreC = r2_score(yVal,yCheckC)
# rScoreS = r2_score(yVal,yCheckS)

# i = 0

# plt.figure(10)
# plt.plot(xP2,FitY2, label = f'Sellmerier fit R^2 = {rScoreS}',linestyle='--')
# plt.plot(xP2,yFit2extrap, label = f'Cauchy fit R^2 = {rScoreC}', linestyle='--')
# plt.scatter(ListExtract(vars()[files[i]+'_Data'],0),ListExtract(vars()[files[i]+'_Data'],13),label = f'Data for {files[i]}', linestyle='--')
# plt.legend()

       
def HiddenCodeForID():
    # plt.figure(0)
    # plt.minorticks_on()
    # plt.grid(which='major', color='k', linestyle='-')
    # plt.grid(which='minor', color='darkgray', linestyle='--')
    # for i in range(len(files)):
    #     plt.errorbar(ListExtract(vars()[files[i]+'_Data'],0),ListExtract(vars()[files[i]+'_Data'],5),ListExtract(vars()[files[i]+'_Data'],6),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1, label = f'n = {np.around(vars()[f"interceptn1-{files[i]}"],3)} + {np.around(vars()[f"mn1-{files[i]}"],3):.2e}/λ^2')
    #     plt.plot(xP,vars()[f"yFitn1 {files[i]}"], '--' , color = LineColour[i])
    # #LabList[i] + f', n =
    # plt.title('Refractive Index of Ga₂O₃ phases (n1)')
    # plt.ylabel('n')
    # plt.xlabel('λ (nm)')
    # plt.legend(fontsize = '10')
    # plt.savefig('n1.png', dpi = 600)
    
    # plt.figure(1)
    # plt.minorticks_on()
    # plt.grid(which='major', color='k', linestyle='-')
    # plt.grid(which='minor', color='darkgray', linestyle='--')
    # for i in range(len(files)):
    #     plt.errorbar(ListExtract(vars()[files[i]+'_Data'],0),ListExtract(vars()[files[i]+'_Data'],13),ListExtract(vars()[files[i]+'_Data'],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1, label = f'n = {np.around(vars()[f"interceptn2-{files[i]}"],3)} + {np.around(vars()[f"mn2-{files[i]}"],3):.2e}/λ^2')
    #     plt.plot(xP,vars()[f"yFitn2 {files[i]}"], '--' , color = LineColour[i] )
    # plt.title('Refractive Index of Ga₂O₃ phases (n2)')
    # plt.ylabel('n')
    # plt.xlabel('λ (nm)')
    # plt.legend(fontsize = '10')
    # plt.savefig('n2.png', dpi = 600)   
    
    
    # num = len(files)
    # numList = list(range(len(files)))
    
    # SubRow,SubCol = choose_subplot_dimensions(num)
    
    # fig, axs = plt.subplots(nrows=SubRow, ncols=SubCol, figsize=(15, 12))
    # plt.subplots_adjust(hspace=0.5)
    # fig.suptitle(f"Cauchy Fit for n1", fontsize=18, y=0.95)
    
    # # loop through tickers and axes
    # for i, ax in zip(numList, axs.ravel()):
    
    #     ax.errorbar(ListExtract(vars()[files[i]+'_Data'],0),ListExtract(vars()[files[i]+'_Data'],13),ListExtract(vars()[files[i]+'_Data'],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1, label = f'n = {np.around(vars()[f"interceptn2-{files[i]}"],3)} + {np.around(vars()[f"mn2-{files[i]}"],3):.2e}/λ^2')
    #     ax.plot(xP,vars()[f"yFitn2 {files[i]}"], '--' , color = LineColour[i] )
    
    # plt.show()
    
    
    
    #loop through the length of tickers and keep track of index
    # num = len(files)
    # plt.subplots_adjust(hspace=1,wspace=1)
    # SubRow,SubCol = choose_subplot_dimensions(num)
    # plt.subplots(figsize=(15,12),sharex=True)
    # plt.suptitle('n1 κ-Ga₂O₃')
    # for i in range(len(files)):
    
    #     # add a new subplot iteratively
    #     plt.minorticks_on()
    #     plt.grid(which='major', color='k', linestyle='-')
    #     plt.grid(which='minor', color='darkgray', linestyle='--')
    #     ax = plt.subplot(SubRow, SubCol, i + 1)
    
    #     # filter df and plot ticker on the new subplot axis
    #     ax.errorbar(ListExtract(vars()[files[i]+'_Data'],0),ListExtract(vars()[files[i]+'_Data'],5),ListExtract(vars()[files[i]+'_Data'],6),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1, label = f'n = {np.around(vars()[f"interceptn1-{files[i]}"],3)} + {np.around(vars()[f"mn1-{files[i]}"],3):.2e}/λ^2')
    #     ax.plot(xP,vars()[f"yFitn1 {files[i]}"], '--' , color = LineColour[i] )
    #     ax.legend()
    #     ax.title.set_text(f'n1 for {files[i]}')
    
    # plt.minorticks_on()
    # plt.grid(which='major', color='k', linestyle='-')
    # plt.grid(which='minor', color='darkgray', linestyle='--')
    # plt.savefig('KappaN1.png', dpi = 600)
    
    
    
    
    
    # plt.figure(2)
    # plt.minorticks_on()
    # plt.grid(which='major', color='k', linestyle='-')
    # plt.grid(which='minor', color='darkgray', linestyle='--')   
    # for i in range(len(files)):
    #     plt.plot((ListExtract(vars()[files3[i]+'_Data'],0)),np.array(ListExtract(vars()[files3[i]+'_Data'],1))*100, label = files[i], color = LineColour[i])
    # plt.title('Transmittance of Ga₂O₃ phases')
    # plt.ylabel('T%')
    # plt.xlabel('λ (nm)')
    # plt.legend(fontsize = '10')
    # plt.savefig('T.png', dpi = 600)

    return