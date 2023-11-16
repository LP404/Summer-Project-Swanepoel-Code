import numpy as np
import copy
import TestArrays as TA
import matplotlib.pyplot as plt


xFitMax,xFitMin,yFitMax,yFitMin,x,y = TA.HereAreArrays()


def MakeTheDammArray(leng):
    weird = np.array([])
    for i in range(leng):
        weird = np.append(weird,1/(np.log(i)+1))
        
    return weird
    


yFitDupMax = copy.deepcopy(yFitMax)
yFitDupMin = copy.deepcopy(yFitMin)

multiplier = MakeTheDammArray(601)

for i in range(len(x)):
    
    if yFitDupMax[i] < y[i]:
    
        corrFactorMax = y[i]/yFitDupMax[i]
                
        yFitDupMax[i] = corrFactorMax*yFitDupMax[i]
        
        
    #     for j in range(len(yFitDupMax[:np.where(yFitDupMax == yFitDupMax[i])[0][0]+1])):

    #         yFitDupMax[j] = corrFactorMax * yFitDupMax[j] * multiplier[j]
            
        
    #     for k in range(len(yFitDupMax[np.where(yFitDupMax == yFitDupMax[i])[0][0]+1:])):
                    
    #         yFitDupMax[k] = corrFactorMax * yFitDupMax[k] * multiplier[k]
    else:
    
        pass
            
# for i in range(len(x)):
    
    if yFitDupMin[i] > y[i]:
    
        corrFactorMin = y[i]/yFitDupMin[i]
        
        yFitDupMin[i] = corrFactorMin*yFitDupMin[i]
        
#         for j in range(len(yFitDupMin[:np.where(yFitDupMin == yFitDupMin[i])[0][0]])):
        
#             yFitDupMin[j] = corrFactorMin*yFitDupMin[j]*np.exp(-1*abs(x[i]-x[j]))
        
#         for k in range(len(yFitDupMin[np.where(yFitDupMin == yFitDupMin[i])[0][0]+1:])):
        
#             yFitDupMin[k] = corrFactorMin*yFitDupMin[k]*np.exp(-1*abs(x[i]-x[k]))


    else:
    
        pass
    
# plt.xlim(250,500)
# plt.ylim(0.7,0.85)
plt.plot(x,y)
plt.plot(xFitMax,yFitMax)
plt.plot(xFitMax,yFitDupMax)
plt.plot(xFitMin,yFitDupMin)