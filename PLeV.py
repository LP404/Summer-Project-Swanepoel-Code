import numpy as np
import nPlotterFunctions as nPF



def ConvertFromNmToEv(x,y):
    h = 6.62607015e-34
    c = 299792458
    
    xNew = ((h*c) / (x*1e-9)) * 6.242e18
    yNew = y * (((x*1e-9)**2)/(h*c))
     
    
    return xNew,yNew


AlphaFileNames, AlphaData = nPF.txtGrabber('PLeV.py','\\PLTXT\\PLeV\\a',False)
BetaFileNames, BetaData = nPF.txtGrabber('PLeV.py','\\PLTXT\\PLeV\\b',False)
KappaFileNames, KappaData = nPF.txtGrabber('PLeV.py','\\PLTXT\\PLeV\\k',False)


for i in range(len(AlphaFileNames)):
    vars()[f'{AlphaFileNames[i]} ev'] = ConvertFromNmToEv(AlphaData[i][0],AlphaData[i][1])
    np.savetxt(f'./PLTXT/PLeV/a/{AlphaFileNames[i]} ev.txt', np.c_[vars()[f'{AlphaFileNames[i]} ev'][0],vars()[f'{AlphaFileNames[i]} ev'][1]], delimiter = "," ,fmt='%f')

    
for i in range(len(BetaFileNames)):
    vars()[f'{BetaFileNames[i]} ev'] = ConvertFromNmToEv(BetaData[i][0],BetaData[i][1])
    np.savetxt(f'./PLTXT/PLeV/b/{BetaFileNames[i]} ev.txt', np.c_[vars()[f'{BetaFileNames[i]} ev'][0],vars()[f'{BetaFileNames[i]} ev'][1]], delimiter = "," ,fmt='%f')

    
for i in range(len(KappaFileNames)):
    vars()[f'{KappaFileNames[i]} ev'] = ConvertFromNmToEv(KappaData[i][0],KappaData[i][1])
    np.savetxt(f'./PLTXT/PLeV/k/{KappaFileNames[i]} ev.txt', np.c_[vars()[f'{KappaFileNames[i]} ev'][0],vars()[f'{KappaFileNames[i]} ev'][1]], delimiter = "," ,fmt='%f')
