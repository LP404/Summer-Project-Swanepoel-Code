import bezier
import numpy as np
import matplotlib.pyplot as plt
import Functions as F
import time


def lineFinder(x1,y1,x2,y2):
    
    m = (y2-y1) / (x2-x1)
    
    c = y1 - (m*x1)
    
    return m,c
    

nodes = np.asfortranarray([[0.0, 0.625, 1.0,1.0],[0.0, 0.5  , 0.5,0.25]])
curve = bezier.Curve(nodes, degree=3)
s_vals = np.linspace(0, 1, 1001)

points = np.linspace(0, 1, 1001)

locX = np.array([])
locY = np.array([])
for i in range(len(points)):
    locx = curve.evaluate_hodograph(points[i])[0][0]
    locy = curve.evaluate_hodograph(points[i])[1][0]
    locX = np.append(locX,locx)
    locY = np.append(locY,locy)

out = curve.evaluate_multi(s_vals)



for i in range(len(out[0])):
    vars()[f'Point X {i}'] = out[0][i] + locX[i]
    vars()[f'Point Y {i}'] = out[1][i] + locY[i]
    vars()[f'm{i}'], vars()[f'c{i}']  = lineFinder(out[0][i],out[1][i],vars()[f'Point X {i}'],vars()[f'Point Y {i}'])

    vars()[f'X {i}'] = np.linspace(0,vars()[f'Point X {i}'],1001)
    vars()[f'Y {i}'] = vars()[f'm{i}']*vars()[f'X {i}'] + vars()[f'c{i}']
    
    
    vars()[f'd {i}'] = np.sqrt(vars()[f'Point X {i}']**2 + vars()[f'Point Y {i}']**2)
    vars()[f'Norm X {i}'] = vars()[f'Point X {i}'] / vars()[f'd {i}']
    vars()[f'Norm Y {i}'] = vars()[f'Point Y {i}'] / vars()[f'd {i}']



for i in range(len(out[0])):
    plt.figure(i)
    plt.scatter(vars()[f'Point X {i}'],vars()[f'Point Y {i}'])
    plt.plot(vars()[f'X {i}'],vars()[f'Y {i}'])
    plt.scatter(vars()[f'Norm X {i}'],vars()[f'Norm Y {i}'])
    plt.plot(out[0],out[1])
    plt.scatter(out[0][i],out[1][i])
    plt.show()
    time.sleep(0.005)

# plt.scatter(point[0],point[1])

# d = np.sqrt(point[0]**2 + point[1]**2)
# NormX = point[0] / d
# NormY = point[1] / d

# plt.scatter(NormX,NormY)
