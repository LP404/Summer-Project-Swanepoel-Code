import numpy as np
import Functions as F

#Factor in substrate
#Add in reflection?

#Mins


m = np.array([])

lambda1 = np.array([100,200,300,400])
n = np.array([2.0,1.99,1.98,1.97])
d = np.array([300,400,500])

LoopVal = len(d) - 1

for i in range(len(lambda1)):
    if i == 0:
        mCalc =  2*n[0]*d[0] / lambda1[0]
        m = np.append(m,np.around(mCalc,0))
    elif i == (len(lambda1)-1):
        mCalc =  2*n[i]*d[i-1] / lambda1[i]
        m = np.append(m,np.around(mCalc,0))
    else:
        for k  in range(0,len(d) - 1):
            mCalc =  2*n[i]*d[i-1+k] / lambda1[i]
            print(d[i-1+k])
            m = np.append(m,np.around(mCalc,0))
        pass
        

newLambda = np.array([])

for i in range(len(lambda1)):
    if i == 0:
        lamCalc =  2*n[0]*d[0] / m[0]
        newLambda = np.append(newLambda,mCalc)
    elif i == (len(lambda1)-1):
        lamCalc =  2*n[i]*d[i-1] / m[i+2]
        newLambda = np.append(newLambda,mCalc)
    else:
        for k  in range(0,len(d) - 1):
            lamCalc =  2*n[i]*d[i-1+k] / m[i+k]
            newLambda = np.append(newLambda,mCalc)
