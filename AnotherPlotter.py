import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import csv
from scipy import stats
from scipy.optimize import curve_fit
import math
import Functions as F
from sklearn.metrics import r2_score
import nPlotterFunctions as nPF
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from decimal import Decimal
import random
import matplotlib as mpl


#Works for positive exponents only
def sci_notation(number, sig_fig=2):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    a, b = ret_string.split("e")
    # remove leading "+" and strip leading zeros
    b = str(int(b))
    return a,b


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


def TickFinder(val):
    if val >= 10:
        return int(np.log10(val))
    else:
        var = float(Decimal(val) % 1)
        if str(val)[::-1].find('.') < str(var)[::-1].find('.'):
            return np.around(var,(str(val)[::-1].find('.')+1))
        else:
            return var

def underline_annotation(text):
    f = plt.gcf()
    ax = plt.gca()
    tb = text.get_tightbbox(f.canvas.get_renderer()).transformed(f.transFigure.inverted())
                            # text isn't drawn immediately and must be 
                            # given a renderer if one isn't cached.
                                                    # tightbbox return units are in 
                                                    # 'figure pixels', transformed 
                                                    # to 'figure fraction'.
    ax.annotate('', xy=(tb.xmin,tb.y0), xytext=(tb.xmax,tb.y0),
                xycoords="figure fraction",
                arrowprops=dict(arrowstyle="-", color='k'))
    
    plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))

axesLabelSize = 32
axesLabelPad = 32
axesFontSize = 26 
axesFontPadding = 20
lineWidth = 5
legendFontSize = 24

h = 6.63e-34
c = 3e8

CLaxesLabelSize = 40
CLaxesLabelPad = 32
CLaxesFontSize = 40
CLaxesFontPadding = 30
CLlineWidth = 5
CLsampleFontSize = 40
CLlegendFontSize = 32

GraphResolution = 100
sizeFigSqaure = CTI(30)

DeltaLam = 0.3

tickMajorSize = 15
tickMinorSize = 10


plt.rcParams["figure.figsize"] = (CTI(15),CTI(10))



h = 6.62607015e-34
c = 299792458


DeltaLam = 0.3


thickness = np.array([350e-7,750e-7,740e-7])
thickness_err = np.array([2.9e-7,2.6e-7,4.5e-7])



path, dirs, files, data, header = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\α','.csv')
path1, dirs1, files1, data1, header1 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\β','.csv')
path2, dirs2, files2, data2, header2 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\κ','.csv')

pathA, dirsA, filesA, dataA, headerA = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\α','.csv')
path1A, dirs1A, files1A, data1A, header1A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\β','.csv')
path2A, dirsA, files2A, data2A, header2A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\κ','.csv')

pathT, dirsT, filesT, dataT, headerT = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\α','.csv')
path1T, dirs1T, files1T, data1T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\β','.csv')
path2T, dirs2T, files2T, data2T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\κ','.csv')





#Ind 0,5,6 n1
#Ind 0,13,14


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple','red','black','gold','orange','salmon','silver','pink','teal','magenta','darkgrey']
ColourComp = ['midnightblue','darkgreen','rebeccapurple']

xP = np.linspace(190,800, 10001)

if data != 0:
    for i in range(len(files)):
        vars()[f"n1yFit {files[i]}"], vars()[f"n2yFit {files[i]}"], vars()[f"n1Cauch {files[i]}"], vars()[f"n2Cauch {files[i]}"], vars()[f"n1CoefList {files[i]}"], vars()[f"n2CoefList {files[i]}"] = nPF.CauchyFinder(files[i],data[i],xP,5,13,0)   
        vars()[f'Lam {filesA[i]}'] = np.array(ListExtract(dataA[i],0))
        vars()[f'Trans {filesA[i]}'] = np.array(ListExtract(dataA[i],1))
        vars()[f'GeoMeanTrans {filesA[i]}'] = np.array(ListExtract(dataA[i],4))
        vars()[f'Alpha {filesA[i]}'] = np.array(ListExtract(dataA[i],9))

else:
    pass

if data1 != 0:
    for i in range(len(files1)):
        vars()[f"n1yFit {files1[i]}"], vars()[f"n2yFit {files1[i]}"], vars()[f"n1Cauch {files1[i]}"], vars()[f"n2Cauch {files1[i]}"], vars()[f"n1CoefList {files1[i]}"], vars()[f"n2CoefList {files1[i]}"] = nPF.CauchyFinder(files1[i],data1[i],xP,5,13,0)
        vars()[f'Lam {files1A[i]}'] = np.array(ListExtract(data1A[i],0))
        vars()[f'Trans {files1A[i]}'] = np.array(ListExtract(data1A[i],1))
        vars()[f'GeoMeanTrans {files1A[i]}'] = np.array(ListExtract(data1A[i],4))
        vars()[f'Alpha {files1A[i]}'] = np.array(ListExtract(data1A[i],9))

else:
    pass

if data2 != 0:
    for i in range(len(files2)):
        vars()[f"n1yFit {files2[i]}"], vars()[f"n2yFit {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2Cauch {files2[i]}"], vars()[f"n2CoefList {files2[i]}"], vars()[f"n2CoefList {files2[i]}"] = nPF.CauchyFinder(files2[i],data2[i],xP,5,13,0)
        vars()[f'Lam {files2A[i]}'] = np.array(ListExtract(data2A[i],0))
        vars()[f'Trans {files2A[i]}'] = np.array(ListExtract(data2A[i],1))
        vars()[f'GeoMeanTrans {files2A[i]}'] = np.array(ListExtract(data2A[i],4))
        vars()[f'Alpha {files2A[i]}'] = np.array(ListExtract(data2A[i],9))

else:
    pass


plt.figure(40,figsize=(CTI(16),CTI(16)), dpi=600)
plt.minorticks_on()
plt.grid(axis='both',which = 'major',linestyle=':', color = 'black')
plt.grid(axis='both',which = 'minor',linestyle=':')
plt.ylabel('n', fontsize = 18)
plt.xlabel(r'$\lambda$ (nm)', fontsize = 18)

label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size 

label_size = 16
mpl.rcParams['ytick.labelsize'] = label_size 
for i in range(len(files)):
        plt.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        Val,Exp = sci_notation(np.around(vars()[f"n2CoefList {files[i]}"][0],3),2)
        plt.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],2)} + ' + r'$\frac{' f"{Val}" r'\times 10^{'f"{Exp}" r'}' '}{λ^2}$')
for i in range(len(files1)):
        plt.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        Val1,Exp1 = sci_notation(np.around(vars()[f"n2CoefList {files1[i]}"][0],3),2)
        plt.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],2)} + ' + r'$\frac{' f"{Val1}" r'\times 10^{'f"{Exp1}" r'}' '}{λ^2}$')
for i in range(len(files2)):
        plt.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+2] ,fmt='o',ms = '4',ecolor='black', elinewidth=3, capsize=1)
        Val2,Exp2 = sci_notation(np.around(vars()[f"n2CoefList {files2[i]}"][0],3),2)
        plt.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+2],label = r'$n_{\kappa}$' + f'= {np.around(vars()[f"n2CoefList {files2[i]}"][1],2)} + ' + r'$\frac{' f"{Val2}" r'\times 10^{'f"{Exp2}" r'}' '}{λ^2}$')
plt.legend(fontsize = 16)

expval = 2

hv = ((h * c) / ( vars()[f'Lam {filesA[0]}'] * 1e-9)) * 6.242e18

LamMaxBgA = 230.7
LamMinBgA = 223


LamMaxBg1A = 240
LamMinBg1A = 228

LamMaxBg2A = 250
LamMinBg2A = 239

LamRangeA = 200
LamRange1A = 200
LamRange2A = 200

hvMinA = ((h * c) / (LamMaxBgA * 1e-9)) * 6.242e18
hvMaxA = ((h * c) / (LamMinBgA * 1e-9)) * 6.242e18

hvMin1A = ((h * c) / (LamMaxBg1A * 1e-9)) * 6.242e18
hvMax1A = ((h * c) / (LamMinBg1A * 1e-9)) * 6.242e18

hvMin2A = ((h * c) / (LamMaxBg2A * 1e-9)) * 6.242e18
hvMax2A = ((h * c) / (LamMinBg2A * 1e-9)) * 6.242e18


for i in range(len(filesA)):
    vars()[f'FinalY {filesA[i]}'] = (vars()[f'Alpha {filesA[i]}'])**expval

    
for i in range(len(files1A)):
    vars()[f'FinalY {files1A[i]}'] = (vars()[f'Alpha {files1A[i]}'])**expval

for i in range(len(files2A)):
    vars()[f'FinalY {files2A[i]}'] = (vars()[f'Alpha {files2A[i]}'])**expval

for i in range(len(filesA)):
    #plt.figure(300 + i)
    # plt.xlim([1.4,3.1])
    # plt.xticks(np.arange(1.4,3.2,0.1))
    # plt.rc('xtick', labelsize=8)
 
  #  Point = np.where(np.gradient(vars()[f'FinalY {filesA[i]}']) == min(np.gradient(vars()[f'FinalY {filesA[i]}'])))[0][0]
    
    m,intercept,err, c_err = F.LineValFinder(hv,vars()[f'FinalY {filesA[i]}'],LamRangeA,'bg',hvMinA,hvMaxA,True)
    yFit = m * hv + intercept
    Contstraint = (yFit >= 0) & (yFit <= max(vars()[f'FinalY {filesA[i]}']))
    line1start = (min(hv),0)
    line1end = (max(hv),0)    
    line2start = (min(hv[Contstraint]),min(yFit[Contstraint]))
    line2end = (max(hv[Contstraint]),max(yFit[Contstraint]))
    vars()[filesA[i]+'Intercept'] = F.LineIntersection((line1start,line1end),(line2start,line2end))   
    vars()[filesA[i]+'hvConst'] = hv[Contstraint]
    vars()[filesA[i]+'yFitConst'] = yFit[Contstraint]
    # plt.plot(hv[Contstraint],yFit[Contstraint],label = 'Intercept = '+str(np.around(vars()[filesA[i]+'Intercept'][0],2)))
    # plt.ylabel('Transmittance, E_bg =' + str(vars()[filesA[i]+'Intercept'])+'eV')        
    vars()[filesA[i]+'hvCOnt'] = hv[Contstraint]
    vars()[filesA[i]+'YfitCon'] = yFit[Contstraint]
    
    PerThickUncert = (thickness_err[0] / thickness[0]) * 100
    LamUncert = (h*c) / (vars()[filesA[i]+'Intercept'][0] * 1.6022e-19) * 1e9
    PerLamUncert = (DeltaLam / LamUncert) * 100
    FinalPerUncert = (2*PerThickUncert) + PerLamUncert
    
    FinalUncert = (FinalPerUncert / 100) * vars()[filesA[i]+'Intercept'][0]
    
    
for i in range(len(files1A)):

    plt.figure(23456)
    m1,intercept1,err1, c_err1 = F.LineValFinder(hv,vars()[f'FinalY {files1A[i]}'],LamRange1A,'bg',hvMin1A,hvMax1A,True)
    yFit1 = m1 * hv + intercept1
    Contstraint1 = (yFit1 >= 0) & (yFit1 <= max(vars()[f'FinalY {files1A[i]}']))
    line1start1 = (min(hv),0)
    line1end1 = (max(hv),0)    
    line2start1 = (min(hv[Contstraint1]),min(yFit1[Contstraint1]))
    line2end1 = (max(hv[Contstraint1]),max(yFit1[Contstraint1]))
    vars()[files1A[i]+'Intercept'] = F.LineIntersection((line1start1,line1end1),(line2start1,line2end1))   
    vars()[files1A[i]+'hvConst'] = hv[Contstraint1]
    vars()[files1A[i]+'yFitConst'] = yFit1[Contstraint1]
    # plt.plot(hv[Contstraint],yFit1[Contstraint],label = 'Intercept = '+str(np.around(vars()[files1A[i]+'Intercept'][0],2)))
    # plt.ylabel('Transmittance, E_bg =' + str(vars()[files1A[i]+'Intercept'])+'eV')        
    vars()[files1A[i]+'hvCOnt'] = hv[Contstraint1]
    vars()[files1A[i]+'YfitCon'] = yFit1[Contstraint1]
    
    PerThickUncert1 = (thickness_err[1] / thickness[1]) * 100
    LamUncert1 = (h*c) / (vars()[files1A[i]+'Intercept'][0] * 1.6022e-19) * 1e9
    PerLamUncert1 = (DeltaLam / LamUncert1) * 100
    FinalPerUncert1 = (2*PerThickUncert1) + PerLamUncert1
    
    FinalUncert1 = (FinalPerUncert1 / 100) * vars()[files1A[i]+'Intercept'][0]
    
for i in range(len(files2A)):

    plt.figure(23456)
    m2,intercept2,err2, c_err2 = F.LineValFinder(hv,vars()[f'FinalY {files2A[i]}'],LamRange2A,'bg',hvMin2A,hvMax2A,True)
    yFit2 = m2 * hv + intercept2
    Contstraint2 = (yFit2 >= 0) & (yFit2 <= max(vars()[f'FinalY {files2A[i]}']))
    line1start2 = (min(hv),0)
    line1end2 = (max(hv),0)    
    line2start2 = (min(hv[Contstraint2]),min(yFit2[Contstraint2]))
    line2end2 = (max(hv[Contstraint2]),max(yFit2[Contstraint2]))
    vars()[files2A[i]+'Intercept'] = F.LineIntersection((line1start2,line1end2),(line2start2,line2end2))   
    vars()[files2A[i]+'hvConst'] = hv[Contstraint2]
    vars()[files2A[i]+'yFitConst'] = yFit2[Contstraint2]
    # plt.plot(hv[Contstraint],yFit2[Contstraint],label = 'Intercept = '+str(np.around(vars()[files2A[i]+'Intercept'][0],2)))
    # plt.ylabel('Transmittance, E_bg =' + str(vars()[files2A[i]+'Intercept'])+'eV')        
    vars()[files2A[i]+'hvCOnt'] = hv[Contstraint2]
    vars()[files2A[i]+'YfitCon'] = yFit2[Contstraint2]
    
    PerThickUncert2 = (thickness_err[2] / thickness[2]) * 100
    LamUncert2 = (h*c) / (vars()[files2A[i]+'Intercept'][0] * 1.6022e-19) * 1e9
    PerLamUncert2 = (DeltaLam / LamUncert2) * 100
    FinalPerUncert2 = (2*PerThickUncert1) + PerLamUncert2
    
    FinalUncert2 = (FinalPerUncert2 / 100) * vars()[files2A[i]+'Intercept'][0]

for i in range(len(filesA)):
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(2*CTI(24),CTI(24)),dpi = GraphResolution)
    fig.subplots_adjust(hspace=0,wspace=0.3)   
    # fig.suptitle('Optical properties of Ga₂O₃ polymorphs', fontsize = 14)
    # fig.minorticks_on()
    # fig.grid(which='major', color='k', linestyle='-')
    # fig.grid(which='minor', color='darkgray', linestyle='--')
    # fig.subplots_adjust(hspace=1,wspace=1)
    # ax1.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=axesFontSize+4)
    # ax1.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=axesFontSize+4)
    # ax1.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= True,labelsize=axesFontSize+4)
    # ax1.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= True,labelsize=axesFontSize+4)
    # ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Alpha {filesA[i]}'])) / 8))
    # ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Alpha {filesA[i]}'])) / 32))
    # ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(hv)),1) / 0.25))
    # ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(hv)),1) / 1))

    # ax2.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top= False,labelsize=axesFontSize-4)
    # ax2.xaxis.set_tick_params(which='minor', size=3.5, width=2, direction='in', top= False,labelsize=axesFontSize-4)
    # ax2.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right= False,labelsize=axesFontSize-4)
    # ax2.yaxis.set_tick_params(which='minor', size=3.5, width=2, direction='in', right= False,labelsize=axesFontSize-4)
    # ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[filesA[i]+'yFitConst'])) / 4))
    # ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[filesA[i]+'yFitConst'])) / 20))
    # ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[filesA[i]+'hvConst'])),1) / 2))
    # ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[filesA[i]+'hvConst'])),1) / 10))
    
    
    MarkSize = 7
    ErrLineWidth = 3
    ErrLineCap = 5

    LabSize = 20
    LegSize = 30
    TickSize = 20
    
    ax1.minorticks_on()
    ax1.grid(which='major', color='k', linestyle='-')
    ax1.grid(which='minor', color='darkgray', linestyle='--')  
    ax1.plot(hv,vars()[f'FinalY {filesA[i]}'], color = LineColour[0])
#    ax1.plot(vars()[filesA[i]+'hvConst'],vars()[filesA[i]+'yFitConst'], label = r'$E_{bg-\alpha}$ = '+str(np.around(vars()[filesA[i]+'Intercept'][0],2))+ r' $\pm$ ' + str(np.around(FinalUncert,2)) + r'eV', linestyle = '--', color = 'xkcd:darkblue')
    ax1.plot(vars()[filesA[i]+'hvConst'],vars()[filesA[i]+'yFitConst'], '--' , lw = 2,label = r'$E_{g-\alpha}$' ,color = ColourComp[0])
    ax1.plot(hv,vars()[f'FinalY {files1A[i]}'], color = LineColour[1])
#    ax1.plot(vars()[files1A[i]+'hvConst'],vars()[files1A[i]+'yFitConst'], label = r'$E_{bg-\beta}$ = '+str(np.around(vars()[files1A[i]+'Intercept'][0],2))+ r' $\pm$ ' + str(np.around(FinalUncert1,2)) + r'eV', linestyle = '--', color = 'xkcd:darkgreen')
    ax1.plot(vars()[files1A[i]+'hvConst'],vars()[files1A[i]+'yFitConst'], '--' , lw = 2,label = r'$E_{g-\beta}$' ,color = ColourComp[1])
    ax1.plot(hv,vars()[f'FinalY {files2A[i]}'], color = LineColour[2])
 #   ax1.plot(vars()[files2A[i]+'hvConst'],vars()[files2A[i]+'yFitConst'], label = r'$E_{bg-\kappa}$ = '+str(np.around(vars()[files2A[i]+'Intercept'][0],2))+ r' $\pm$ ' + str(np.around(FinalUncert2,2)) + r'eV', linestyle = '--', color = 'xkcd:plum')
    ax1.plot(vars()[files2A[i]+'hvConst'],vars()[files2A[i]+'yFitConst'], '--' , lw = 2,label = r'$E_{g-\kappa}$' ,color = ColourComp[2])
    ax1.set_xlim(4.5,5.6)
    ax1.set_ylim(0,vars()[filesA[i]+'yFitConst'][np.where(np.around(vars()[filesA[i]+'hvConst'],2) == 5.45)[0][0]])
    ax1.set_ylabel(r'$\alpha^{2} \, (cm^{-1})^{2}$ (A.U.)', fontsize = LabSize)
    ax1.set_xlabel(r'h$\nu$ (eV)', fontsize = LabSize)
    ax1.tick_params(axis='both', which='major', labelsize=TickSize)
    ax1.legend(fontsize = LegSize)



    ax2.minorticks_on()
    ax2.grid(which='major', color='k', linestyle='-')
    ax2.grid(which='minor', color='darkgray', linestyle='--')
    ax2.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = MarkSize,ecolor='black', elinewidth=ErrLineWidth, capsize=ErrLineCap)
#    ax2.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],2)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/λ^2')
    ax2.plot(xP, vars()[f"n2yFit {files[i]}"], '--',color = LineColour[i],label = r'$n_{\alpha}$')
    ax2.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = MarkSize,ecolor='black', elinewidth=ErrLineWidth, capsize=ErrLineCap)
#    ax2.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],2)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/λ^2')
    ax2.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',color = LineColour[i+1],label = r'$n_{\beta}$')
    ax2.errorbar(ListExtract(data2[i],0),ListExtract(data2[i],13),ListExtract(data2[i],14),DeltaLam, color = LineColour[i+2] ,fmt='o',ms = MarkSize,ecolor='black', elinewidth=ErrLineWidth, capsize=ErrLineCap)
#    ax2.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+2],label = r'$n_{\kappa}$' + f'= {np.around(vars()[f"n2CoefList {files2[i]}"][1],2)} + {np.around(vars()[f"n2CoefList {files2[i]}"][0],3):.2e}/λ^2')
    ax2.plot(xP, vars()[f"n2yFit {files2[i]}"], '--',color = LineColour[i+2],label = r'$n_{\kappa}$')
    ax2.set_ylabel(r'Refractive Index', fontsize = LabSize)
    ax2.set_xlabel(r'$\lambda$ (nm)', fontsize = LabSize)
    ax2.tick_params(axis='both', which='major', labelsize=TickSize)
    ax2.legend(fontsize = LegSize)


