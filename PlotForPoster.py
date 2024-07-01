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
axesLabelPad = 10
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

path, dirs, files, data, header = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\α','.csv')
path1, dirs1, files1, data1, header1 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\β','.csv')
path2, dirs2, files2, data2, header2 = nPF.fileGrabberTrans('BandGapCalc.py','\\CSVout\\Transmittance\\κ','.csv')

pathA, dirsA, filesA, dataA, headerA = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\α','.csv')
path1A, dirs1A, files1A, data1A, header1A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\β','.csv')
path2A, dirsA, files2A, data2A, header2A = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\κ','.csv')


pathSub, dirsSub, filesSub, dataSub, headerSub = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Absorption\\Sub','.csv')


pathT, dirsT, filesT, dataT, headerT = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\α','.csv')
path1T, dirs1T, files1T, data1T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\β','.csv')
path2T, dirs2T, files2T, data2T, header2T = nPF.fileGrabber('BandGapCalc.py','\\CSVout\\Thickness\\κ','.csv')


headerAlphaCL, dataAlphaCL = nPF.datGrabber('ExportedDataAlpha','fitykExports\\')
headerBetaCL, dataBetaCL = nPF.datGrabber('ExportedDataBeta','fitykExports\\')
headerKappaCL, dataKappaCL = nPF.datGrabber('ExportedDataKappa','fitykExports\\')

headerAlphaCLLT, dataAlphaCLLT = nPF.datGrabber('MeanSpectra10keV','fitykExports\\')
TheaderBetaCLLT, dataBetaCLLT = nPF.datGrabber('BetaMeanSpectra','fitykExports\\')
headerKappaCLLT, dataKappaCLLT = nPF.datGrabber('MeanKappaSpecrtra','fitykExports\\')

suffix = ".csv"

#Ind 0,5,6 n1
#Ind 0,13,14


LabList = ['α-Ga₂O₃','β-Ga₂O₃','κ-Ga₂O₃']
LineColour = ['blue','green','purple','red','black','gold','orange','salmon','silver','pink','teal','magenta','darkgrey']
LineStyle = ['solid','dashed']

xP = np.linspace(190,800, 10001)



if dataSub != 0:
    for i in range(len(filesSub)):
        vars()[f'Lam {filesSub[i]}'] = np.array(ListExtract(dataSub[i],0))
        vars()[f'Trans {filesSub[i]}'] = np.array(ListExtract(dataSub[i],1))

else:
    pass



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


# fig,ax = plt.subplots(figsize = (sizeFigSqaure,sizeFigSqaure),dpi = GraphResolution)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(True)
# ax.spines['bottom'].set_visible(True)
# ax.spines['top'].set_visible(False)
# ax.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=axesFontSize)
# ax.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=axesFontSize)
# ax.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=axesFontSize)
# ax.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=axesFontSize)
# ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(xP))))
# ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(xP)) / 5))
# ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f"n2yFit {files[0]}"])),1) / 5))
# ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f"n2yFit {files[0]}"])),1) / 25))
# #ax.set_title('n values of Ga₂O₃ polymorphs', pad=10)
# ax.text(0.85, 0.29, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, color = LineColour[0], fontsize = legendFontSize)
# ax.text(0.85, 0.11, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, color = LineColour[1], fontsize = legendFontSize)
# ax.set_ylabel(r'n', labelpad = axesLabelPad, fontsize = axesLabelSize)
# ax.set_xlabel(r'$\lambda$ (nm)', labelpad = axesLabelPad, fontsize = axesLabelSize)
# for i in range(len(files)):
#         ax.errorbar(ListExtract(data[i],0),ListExtract(data[i],13),ListExtract(data[i],14),DeltaLam, color = LineColour[i] ,fmt='o',ms = '12',ecolor='black', elinewidth=8, capsize=1)
#         ax.plot(xP, vars()[f"n2yFit {files[i]}"], '--',linewidth = lineWidth,color = LineColour[i],label = r'$n_{\alpha}$' + f'= {np.around(vars()[f"n2CoefList {files[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files[i]}"][0],3):.2e}/$λ^{2}$')
# for i in range(len(files1)):
#         ax.errorbar(ListExtract(data1[i],0),ListExtract(data1[i],13),ListExtract(data1[i],14),DeltaLam, color = LineColour[i+1] ,fmt='o',ms = '12',ecolor='black', elinewidth=8, capsize=1)
#         ax.plot(xP, vars()[f"n2yFit {files1[i]}"], '--',linewidth = lineWidth,color = LineColour[i+1],label = r'$n_{\beta}$' + f'= {np.around(vars()[f"n2CoefList {files1[i]}"][1],3)} + {np.around(vars()[f"n2CoefList {files1[i]}"][0],3):.2e}/$λ^{2}$')
# ax.tick_params(axis='both', which='major', pad=axesFontPadding)
# ax.legend(bbox_to_anchor=(0.95, 0.95), fontsize=legendFontSize)
# plt.savefig("RI Poster.png",bbox_inches="tight")

Thickness = np.array([352,754,852])

hv = ((h * c) / ( vars()[f'Lam {filesA[0]}'] * 1e-9)) * 6.242e18

### Values for exp of 2

expval = 2

LamMaxBgA = 230.7
LamMinBgA = 223

LamMaxBg1A = 240
LamMinBg1A = 228

LamMaxBg2A = 246
LamMinBg2A = 240

LamRangeA = 200
LamRange1A = 200
LamRange2A = 200

hvMinA = ((h * c) / (LamMaxBgA * 1e-9)) * 6.242e18
hvMaxA = ((h * c) / (LamMinBgA * 1e-9)) * 6.242e18

hvMin1A = ((h * c) / (LamMaxBg1A * 1e-9)) * 6.242e18
hvMax1A = ((h * c) / (LamMinBg1A * 1e-9)) * 6.242e18

hvMin2A = ((h * c) / (LamMaxBg2A * 1e-9)) * 6.242e18
hvMax2A = ((h * c) / (LamMinBg2A * 1e-9)) * 6.242e18



### Values for exp of 0.5
###Dont't work

# expval = 0.5

# LamMaxBgA = 234
# LamMinBgA = 222

# LamMaxBg1A = 260
# LamMinBg1A = 245

# LamMaxBg2A = 251
# LamMinBg2A = 245

# LamRangeA = 200
# LamRange1A = 200
# LamRange2A = 200

# hvMinA = ((h * c) / (LamMaxBgA * 1e-9)) * 6.242e18
# hvMaxA = ((h * c) / (LamMinBgA * 1e-9)) * 6.242e18

# hvMin1A = ((h * c) / (LamMaxBg1A * 1e-9)) * 6.242e18
# hvMax1A = ((h * c) / (LamMinBg1A * 1e-9)) * 6.242e18

# hvMin2A = ((h * c) / (LamMaxBg2A * 1e-9)) * 6.242e18
# hvMax2A = ((h * c) / (LamMinBg2A * 1e-9)) * 6.242e18


# for i in range(len(filesA)):
#     vars()[f'FinalY {filesA[i]}'] = (vars()[f'Alpha {filesA[i]}'])**expval

    
# for i in range(len(files1A)):
#     vars()[f'FinalY {files1A[i]}'] = (vars()[f'Alpha {files1A[i]}'])**expval

# for i in range(len(files1A)):
#     vars()[f'FinalY {files2A[i]}'] = (vars()[f'Alpha {files2A[i]}'])**expval

    
for i in range(len(filesA)):
    
    vars()[f'FinalT {filesA[i]}'] =   vars()[f'GeoMeanTrans {filesA[i]}'] / vars()[f'Trans {filesSub[i]}']
    
    
    vars()[f'Final Alpha {filesA[i]}'] = (np.log((1/vars()[f'FinalT {filesA[i]}'])) / Thickness[0] ) * 1e7
    
    vars()[f'FinalY {filesA[i]}'] = (vars()[f'Final Alpha {filesA[i]}'])**expval

    
for i in range(len(files1A)):
    
    vars()[f'FinalT {files1A[i]}'] =   vars()[f'GeoMeanTrans {files1A[i]}'] / vars()[f'Trans {filesSub[i]}']
    
    
    vars()[f'Final Alpha {files1A[i]}'] = (np.log((1/vars()[f'FinalT {files1A[i]}'])) / Thickness[1]) * 1e7
    
    vars()[f'FinalY {files1A[i]}'] = (vars()[f'Final Alpha {files1A[i]}'])**expval

for i in range(len(files2A)):
    
    vars()[f'FinalT {files2A[i]}'] =   vars()[f'GeoMeanTrans {files2A[i]}'] / vars()[f'Trans {filesSub[i]}']
    
    
    vars()[f'Final Alpha {files2A[i]}'] = (np.log((1/vars()[f'FinalT {files2A[i]}'])) / Thickness[2]) * 1e7
    
    vars()[f'FinalY {files2A[i]}'] = (vars()[f'Final Alpha {files2A[i]}'])**expval

for i in range(len(filesA)):
    #plt.figure(300 + i)
    # plt.xlim([1.4,3.1])
    # plt.xticks(np.arange(1.4,3.2,0.1))
    # plt.rc('xtick', labelsize=8)
 
  #  Point = np.where(np.gradient(vars()[f'FinalY {filesA[i]}']) == min(np.gradient(vars()[f'FinalY {filesA[i]}'])))[0][0]
    
    m,intercept,err,intererr = F.LineValFinder(hv,vars()[f'FinalY {filesA[i]}'],LamRangeA,'bg',hvMinA,hvMaxA,True)
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

    
for i in range(len(files1A)):

    plt.figure(23456)
    m1,intercept1,err1,intererr1 = F.LineValFinder(hv,vars()[f'FinalY {files1A[i]}'],LamRange1A,'bg',hvMin1A,hvMax1A,True)
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
 
for i in range(len(files2A)):

    plt.figure(23457)
    m2,intercept2,err2,intererr2 = F.LineValFinder(hv,vars()[f'FinalY {files2A[i]}'],LamRange2A,'bg',hvMin2A,hvMax2A,True)
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


ColourComp = ['midnightblue','darkgreen','rebeccapurple']

for i in range(len(filesA)):
    fig1,ax1 = plt.subplots(figsize = (sizeFigSqaure,sizeFigSqaure),dpi = GraphResolution)
    ax1.spines['right'].set_visible(True)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['top'].set_visible(False)
    ax1.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=axesFontSize+4)
    ax1.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=axesFontSize+4)
    ax1.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= True,labelsize=axesFontSize+4)
    ax1.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= True,labelsize=axesFontSize+4)
    ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Alpha {filesA[i]}'])) / 8))
    ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Alpha {filesA[i]}'])) / 32))
    ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(hv)),1) / 0.25))
    ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(hv)),1) / 1))
    
    ax2 = fig1.add_axes([0.13,0.48,0.4,0.4])
    ax2.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top= False,labelsize=axesFontSize-4)
    ax2.xaxis.set_tick_params(which='minor', size=3.5, width=2, direction='in', top= False,labelsize=axesFontSize-4)
    ax2.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right= False,labelsize=axesFontSize-4)
    ax2.yaxis.set_tick_params(which='minor', size=3.5, width=2, direction='in', right= False,labelsize=axesFontSize-4)
    ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[filesA[i]+'yFitConst'])) / 4))
    ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[filesA[i]+'yFitConst'])) / 20))
    ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[filesA[i]+'hvConst'])),1) / 2))
    ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[filesA[i]+'hvConst'])),1) / 10))

    
    ax1.semilogy(hv,vars()[f'Alpha {filesA[i]}'], label = LabList[0], color = LineColour[0], linewidth = lineWidth)
    ax1.semilogy(hv,vars()[f'Alpha {files1A[i]}'], label = LabList[1], color = LineColour[1], linewidth = lineWidth)
    ax1.semilogy(hv,vars()[f'Alpha {files2A[i]}'], label = LabList[2], color = LineColour[2], linewidth = lineWidth)
    ax1.set_xlabel(r'h$\nu$ (eV)', labelpad = axesLabelPad, fontsize = axesLabelSize+4)
    ax1.set_ylabel(r'$\alpha \, (cm^{-1})$', labelpad = axesLabelPad, fontsize = axesLabelSize+4)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.text(1.15 / (sizeFigSqaure/8), 1.42 / (sizeFigSqaure/8), r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes, color = LineColour[0], fontsize = legendFontSize)
    ax1.text(1.10/ (sizeFigSqaure/8), 0.3 / (sizeFigSqaure/8), r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes, color = LineColour[1], fontsize = legendFontSize)
    ax1.text(0.75/ (sizeFigSqaure/8), 0.35 / (sizeFigSqaure/8), r'$\kappa-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes, color = LineColour[2], fontsize = legendFontSize)

    ax1.tick_params(axis='both', which='major', pad=axesFontPadding)
  

    ax2.plot(hv,vars()[f'FinalY {filesA[i]}'], color = LineColour[0])
    ax2.plot(vars()[filesA[i]+'hvConst'],vars()[filesA[i]+'yFitConst'], label = r'$E_{bg}$ = '+str(np.around(vars()[filesA[i]+'Intercept'][0],2))+r'eV', linestyle = '--', color = ColourComp[0])
    ax2.plot(hv,vars()[f'FinalY {files1A[i]}'], color = LineColour[1])
    ax2.plot(vars()[files1A[i]+'hvConst'],vars()[files1A[i]+'yFitConst'], label = r'$E_{bg}$ = '+str(np.around(vars()[files1A[i]+'Intercept'][0],2))+r'eV', linestyle = '--', color = ColourComp[1])
    ax2.plot(hv,vars()[f'FinalY {files2A[i]}'], color = LineColour[2])
    ax2.plot(vars()[files2A[i]+'hvConst'],vars()[files2A[i]+'yFitConst'], label = r'$E_{bg}$ = '+str(np.around(vars()[files2A[i]+'Intercept'][0],2))+r'eV', linestyle = '--', color = ColourComp[2])
    ax2.set_xlim(4.5,5.6)
    ax2.set_ylim(0,vars()[filesA[i]+'yFitConst'][np.where(np.around(vars()[filesA[i]+'hvConst'],2) == 5.45)[0][0]])
    ax2.set_ylabel(r'$\alpha^{2} \, (cm^{-1})^{2} \times 10^{11}$', labelpad = axesLabelPad-4, fontsize = axesLabelSize-4)
    ax2.set_xlabel(r'h$\nu$ (eV)', labelpad = axesLabelPad-4, fontsize = axesLabelSize-4)
    ax2.legend(bbox_to_anchor=(0.65, 1), fontsize=legendFontSize-4)
    ax2.tick_params(axis='both', which='major', pad=axesFontPadding)
    ax2.yaxis.offsetText.set_color('#FFFFFF')
    plt.savefig("Abs + Ebg Poster.png",bbox_inches="tight")

# plt.figure(20)
# plt.minorticks_on()
# plt.tight_layout()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# plt.semilogy(hv,vars()[f'Alpha {filesA[i]}'], label = LabList[0], color = LineColour[0])
# plt.semilogy(hv,vars()[f'Alpha {files1A[i]}'], label = LabList[1], color = LineColour[1])
# plt.semilogy(hv,vars()[f'Alpha {files2A[i]}'], label = LabList[2], color = LineColour[2])
# plt.xlabel(r'h$\nu$ (eV)', fontsize = 14)
# plt.ylabel(r'$\alpha \, (cm^{-1})$', fontsize = 14)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.legend()

# plt.figure(21)
# plt.minorticks_on()
# plt.tight_layout()
# plt.grid(which='major', color='k', linestyle='-')
# plt.grid(which='minor', color='darkgray', linestyle='--')
# plt.plot(hv,vars()[f'FinalY {filesA[i]}'], color = LineColour[0])
# plt.plot(vars()[filesA[i]+'hvConst'],vars()[filesA[i]+'yFitConst'], label = r'$\alpha - E_{bg}  $', linestyle = '--', color = LineColour[3])
# plt.plot(hv,vars()[f'FinalY {files1A[i]}'], color = LineColour[1])
# plt.plot(vars()[files1A[i]+'hvConst'],vars()[files1A[i]+'yFitConst'], label = r'$\beta - E_{bg} $', linestyle = '--', color = 'darkred')
# plt.plot(hv,vars()[f'FinalY {files2A[i]}'], color = LineColour[2])
# plt.plot(vars()[files2A[i]+'hvConst'],vars()[files2A[i]+'yFitConst'], label = r'$\kappa - E_{bg} $', linestyle = '--', color = 'darkblue')
# plt.xlim(4.5,5.6)
# plt.ylim(0,vars()[filesA[i]+'yFitConst'][np.where(np.around(vars()[filesA[i]+'hvConst'],2) == 5.45)[0][0]])
# plt.ylabel(r'$\alpha^{2} \, (cm^{-1})^{2}$',fontsize = 14)
# plt.xlabel(r'h$\nu$ (eV)', fontsize = 14)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.legend(fontsize = 14)

    fig2,ax4 = plt.subplots(figsize = (sizeFigSqaure,sizeFigSqaure),dpi = GraphResolution)
    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_visible(True)
    ax4.spines['bottom'].set_visible(True)
    ax4.spines['top'].set_visible(False)
    ax4.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=axesFontSize)
    ax4.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=axesFontSize)
    ax4.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=axesFontSize)
    ax4.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=axesFontSize)
    ax4.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Lam {filesA[i]}']))))
    ax4.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(vars()[f'Lam {filesA[i]}'])) / 5))
    ax4.yaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f'GeoMeanTrans {filesA[i]}'])),0) / 10))
    ax4.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(vars()[f'GeoMeanTrans {filesA[i]}'])),0) / 50))
    #ax4.set_title('n values of Ga₂O₃ polymorphs', pad=10)
    # ax4.text(0.85, 0.29, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, color = LineColour[0], fontsize = 14)
    # ax4.text(0.85, 0.11, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, color = LineColour[1], fontsize = 14)
    ax4.set_xlim(190,500)
    ax4.set_ylabel(r'T%', labelpad = axesLabelPad, fontsize = axesLabelSize)
    ax4.set_xlabel(r'$\lambda$ (nm)', labelpad = axesLabelPad, fontsize = axesLabelSize)
    ax4.plot(vars()[f'Lam {filesA[i]}'], vars()[f'GeoMeanTrans {filesA[i]}'],color = LineColour[0], label = 'a', linewidth = lineWidth)
    ax4.plot(vars()[f'Lam {files1A[i]}'], vars()[f'GeoMeanTrans {files1A[i]}'],color = LineColour[1], label = 'b', linewidth = lineWidth)
    ax4.tick_params(axis='both', which='major', pad=axesFontPadding)
    ax4.text(0.8, 0.95, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes, color = LineColour[0], fontsize = legendFontSize)
    ax4.text(0.8,0.875, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes, color = LineColour[1], fontsize = legendFontSize)
    plt.savefig("Transmittance Poster.png",bbox_inches="tight")

# CLpeakColours = ['blue','indigo','violet','black','gold']


# alphaCLpeakNo = headerAlphaCL[-1]
# betaCLpeakNo = headerBetaCL[-1]
# kappaCLpeakNo = headerKappaCL[-1]

# peakEdgeCutoff = 0

# fig = plt.figure(81, figsize = (sizeFigSqaure * 2,sizeFigSqaure),dpi = GraphResolution)
# ax5 = fig.add_axes([0,0,1,1])
# ax5.spines['right'].set_visible(False)
# ax5.spines['left'].set_visible(True)
# ax5.spines['bottom'].set_visible(True)
# ax5.spines['top'].set_visible(False)
# ax5.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax5.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax5.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax5.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax5.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataAlphaCL[5])) / 5))
# ax5.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataAlphaCL[5])) / 25))
# ax5.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataAlphaCL[0])),1) / 2))
# ax5.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataAlphaCL[0])),1) / 8))

# ax5.set_ylabel(r'counts / eV  $\times 10^6$', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax5.set_xlabel(r'$h \nu$ (eV)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax5.set_xlim(1.8,4.6)
# ax5.plot(dataAlphaCL[0][np.where(dataAlphaCL[4] > peakEdgeCutoff)[0]], dataAlphaCL[4][np.where(dataAlphaCL[4] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[0], label = f"Peak at {np.around(dataAlphaCL[0][np.where(dataAlphaCL[4] == max(dataAlphaCL[4]))[0][0]],2)}eV")
# ax5.plot(dataAlphaCL[0][np.where(dataAlphaCL[5] > peakEdgeCutoff)[0]], dataAlphaCL[5][np.where(dataAlphaCL[5] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[1], label = f"Peak at {np.around(dataAlphaCL[0][np.where(dataAlphaCL[5] == max(dataAlphaCL[5]))[0][0]],2)}eV")
# ax5.plot(dataAlphaCL[0][np.where(dataAlphaCL[6] > peakEdgeCutoff)[0]], dataAlphaCL[6][np.where(dataAlphaCL[6] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[2], label = f"Peak at {np.around(dataAlphaCL[0][np.where(dataAlphaCL[6] == max(dataAlphaCL[6]))[0][0]],2)}eV")
# ax5.plot(dataAlphaCL[0], dataAlphaCL[1],color = CLpeakColours[3], linewidth = CLlineWidth, label = "Raw Data")
# ax5.plot(dataAlphaCL[0], dataAlphaCL[7],color = CLpeakColours[4], linewidth = CLlineWidth, label = "Peak Convolution", linestyle = '--')
# ax5.yaxis.offsetText.set_color('#FFFFFF')
# ax5.text(0.7, 0.9, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax5.transAxes, color = LineColour[0], fontsize = CLsampleFontSize)
# plt.legend(fontsize = CLlegendFontSize)
# ax5.set_ylim(bottom = 0)
# ax5.tick_params(axis='both', which='major', pad=CLaxesFontPadding)
# plt.savefig("Alpha CL.png",bbox_inches="tight")


# fig = plt.figure(82, figsize = (sizeFigSqaure * 2 , sizeFigSqaure),dpi = GraphResolution)
# ax6 = fig.add_axes([0,0,1,1])
# ax6.spines['right'].set_visible(False)
# ax6.spines['left'].set_visible(True)
# ax6.spines['bottom'].set_visible(True)
# ax6.spines['top'].set_visible(False)
# ax6.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax6.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax6.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax6.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)

# ax6.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataBetaCL[4])) / 1 ))
# ax6.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataBetaCL[4])) / 5))
# ax6.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataBetaCL[0])),1) / 2))
# ax6.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataBetaCL[0])),1) / 8))


# ax6.set_xlim(1.8,4.2)
# ax6.set_ylabel(r'counts / eV', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax6.set_xlabel(r'$h \nu$ (eV)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax6.plot(dataBetaCL[0][np.where(dataBetaCL[4] > peakEdgeCutoff)[0]], dataBetaCL[4][np.where(dataBetaCL[4] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[1], label = f"Peak at {np.around(dataBetaCL[0][np.where(dataBetaCL[4] == max(dataBetaCL[4]))[0][0]],2)}eV")
# ax6.plot(dataBetaCL[0][np.where(dataBetaCL[5] > peakEdgeCutoff)[0]], dataBetaCL[5][np.where(dataBetaCL[5] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[0], label = f"Peak at {np.around(dataBetaCL[0][np.where(dataBetaCL[5] == max(dataBetaCL[5]))[0][0]],2)}eV")
# ax6.plot(dataBetaCL[0], dataBetaCL[1],color = CLpeakColours[3], linewidth = CLlineWidth, label = "Raw Data")
# ax6.plot(dataBetaCL[0], dataBetaCL[6],color = CLpeakColours[4], linewidth = CLlineWidth, label = "Peak Convolution", linestyle = '--')
# ax6.text(0.7, 0.9, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax6.transAxes, color = LineColour[1], fontsize = CLsampleFontSize)
# ax6.set_ylim(bottom = 0)
# ax6.tick_params(axis='both', which='major', pad=CLaxesFontPadding)
# plt.legend(fontsize = CLlegendFontSize)
# plt.savefig("Beta CL.png",bbox_inches="tight")

# fig = plt.figure(83, figsize = (sizeFigSqaure * 2, sizeFigSqaure),dpi = GraphResolution)
# ax7 = fig.add_axes([0,0,1,1])
# ax7.spines['right'].set_visible(False)
# ax7.spines['left'].set_visible(True)
# ax7.spines['bottom'].set_visible(True)
# ax7.spines['top'].set_visible(False)
# ax7.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax7.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax7.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax7.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)

# ax7.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCL[6])) / 5 ))
# ax7.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCL[6])) / 25))
# ax7.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCL[0])),1) / 2))
# ax7.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCL[0])),1) / 8))

# ax7.set_ylabel(r'counts / eV', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax7.set_xlabel(r'$h \nu$ (eV)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax7.set_xlim(1.8,4)
# ax7.plot(dataKappaCL[0][np.where(dataKappaCL[4] > peakEdgeCutoff)[0]], dataKappaCL[4][np.where(dataKappaCL[4] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[0], label = f"Peak at {np.around(dataKappaCL[0][np.where(dataKappaCL[4] == max(dataKappaCL[4]))[0][0]],2)}eV")
# ax7.plot(dataKappaCL[0][np.where(dataKappaCL[5] > peakEdgeCutoff)[0]], dataKappaCL[5][np.where(dataKappaCL[5] > peakEdgeCutoff)[0]], linewidth = CLlineWidth,color = CLpeakColours[1], label = f"Peak at {np.around(dataKappaCL[0][np.where(dataKappaCL[5] == max(dataKappaCL[5]))[0][0]],2)}eV")
# ax7.plot(dataKappaCL[0], dataKappaCL[1],color = CLpeakColours[3], linewidth = CLlineWidth, label = "Raw Data")
# ax7.plot(dataKappaCL[0], dataKappaCL[6],color = CLpeakColours[4], linewidth = CLlineWidth, label = "Peak Convolution", linestyle = '--')
# ax7.text(0.7, 0.9, r'$\kappa-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax7.transAxes, color = LineColour[2], fontsize = CLsampleFontSize)
# ax7.set_ylim(bottom = 0)
# ax7.tick_params(axis='both', which='major', pad=CLaxesFontPadding)
# plt.legend(fontsize = CLlegendFontSize)
# plt.savefig("Kappa CL.png",bbox_inches="tight")

# fig = plt.figure(83, figsize = (sizeFigSqaure * 2, sizeFigSqaure),dpi = GraphResolution)
# ax8 = fig.add_axes([0,0,1,1])
# ax8.spines['right'].set_visible(False)
# ax8.spines['left'].set_visible(True)
# ax8.spines['bottom'].set_visible(True)
# ax8.spines['top'].set_visible(False)
# ax8.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax8.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax8.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax8.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)

# ax8.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCL[1]/max(dataKappaCL[1]))) / 5 ))
# ax8.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCL[1]/max(dataKappaCL[1]))) / 25))
# ax8.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCL[0])),1) / 2))
# ax8.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCL[0])),1) / 10))

# ax8.set_ylabel(r'CL Intensity (a.u.)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax8.set_xlabel(r'$h \nu$ (eV)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax8.set_xlim(1.8,4.2)
# ax8.plot(dataAlphaCL[0], dataAlphaCL[1]/max(dataAlphaCL[1]),color =  LineColour[0], linewidth = CLlineWidth, linestyle = (0, (3,3)))
# ax8.plot(dataBetaCL[0], dataBetaCL[1]/max(dataBetaCL[1]),color =  LineColour[1], linewidth = CLlineWidth, linestyle = (0,(1,1)))
# ax8.plot(dataKappaCL[0], dataKappaCL[1]/max(dataKappaCL[1]),color =  LineColour[2], linewidth = CLlineWidth)
# leg1 = ax8.text(0.4, 0.9, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax8.transAxes, color = LineColour[0], fontsize = CLsampleFontSize)
# leg2 = ax8.text(0.7, 0.9, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax8.transAxes, color = LineColour[1], fontsize = CLsampleFontSize)
# leg3 = ax8.text(0.15, 0.9, r'$\kappa-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax8.transAxes, color = LineColour[2], fontsize = CLsampleFontSize)
# ax8.set_ylim(bottom = 0)
# ax8.tick_params(axis='both', which='major', pad=CLaxesFontPadding)
# plt.savefig("Spectra CL.png",bbox_inches="tight")

# fig = plt.figure(84, figsize = (sizeFigSqaure * 2, sizeFigSqaure),dpi = GraphResolution)
# ax9 = fig.add_axes([0,0,1,1])
# ax9.spines['right'].set_visible(False)
# ax9.spines['left'].set_visible(True)
# ax9.spines['bottom'].set_visible(True)
# ax9.spines['top'].set_visible(False)
# ax9.xaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax9.xaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', top= False,labelsize=CLaxesFontSize)
# ax9.yaxis.set_tick_params(which='major', size = tickMajorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)
# ax9.yaxis.set_tick_params(which='minor', size = tickMinorSize, width=2, direction='in', right= False,labelsize=CLaxesFontSize)

# ax9.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCLLT[1]/max(dataKappaCLLT[1]))) / 5 ))
# ax9.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10**TickFinder(max(dataKappaCLLT[1]/max(dataKappaCLLT[1]))) / 25))
# ax9.xaxis.set_major_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCLLT[0])),1) / 2))
# ax9.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(np.around(TickFinder(max(dataKappaCLLT[0])),1) / 10))

# ax9.set_ylabel(r'CL Intensity (a.u.)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax9.set_xlabel(r'$h \nu$ (eV)', labelpad = CLaxesLabelPad, fontsize = CLaxesLabelSize)
# ax9.set_xlim(1.8,4.4)
# ax9.plot(dataAlphaCLLT[0], dataAlphaCLLT[1]/max(dataAlphaCLLT[1]),color =  LineColour[0], linewidth = CLlineWidth, linestyle = (0, (3,3)))
# ax9.plot(dataBetaCLLT[0], dataBetaCLLT[1]/max(dataBetaCLLT[1]),color =  LineColour[1], linewidth = CLlineWidth, linestyle = (0,(1,1)))
# ax9.plot(dataKappaCLLT[0], dataKappaCLLT[1]/max(dataKappaCLLT[1]),color =  LineColour[2], linewidth = CLlineWidth)
# leg1 = ax9.text(0.38, 0.9, r'$\alpha-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax9.transAxes, color = LineColour[0], fontsize = CLsampleFontSize)
# leg2 = ax9.text(0.7, 0.9, r'$\beta-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax9.transAxes, color = LineColour[1], fontsize = CLsampleFontSize)
# leg3 = ax9.text(0.18, 0.9, r'$\kappa-Ga_{2}O_{3}$', horizontalalignment='center',verticalalignment='center', transform=ax9.transAxes, color = LineColour[2], fontsize = CLsampleFontSize)
# ax9.set_ylim(bottom = 0)
# ax9.tick_params(axis='both', which='major', pad=CLaxesFontPadding)
# plt.savefig("Spectra LT CL.png",bbox_inches="tight")
