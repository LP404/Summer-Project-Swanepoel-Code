import numpy as np
import Functions as F

#Factor in substrate
#Add in reflection?

#Mins


m = np.array([])

lambda1 = np.array([100,200,300,400])
n = np.array([2.0,1.99,1.98,1.977])
d = np.array([300,400,500])

HasLoop = 0
for i in range(len(lambda1)):
    if i == 0:
        mCalc =  2*n[0]*d[0] / lambda1[0]
        m = np.append(m,np.around(mCalc,0))
    elif i == (len(lambda1)-1):
        mCalc =  2*n[i]*d[i-1] / lambda1[i]
        m = np.append(m,np.around(mCalc,0))
    else:
        for k  in range(0,2):
            mCalc =  2*n[i]*d[i-1+k] / lambda1[i]
            print('Outer Loop = ' + str(i))
            print('Inner Loop = ' + str(k))
            print('Thickness Index Pointer = ' + str(i-1+k))
            print('Guess Index Pointer = ' + str(i+k+HasLoop))
            print(' ')
            m = np.append(m,np.around(mCalc,0))
            if k == 1: 
                HasLoop +=1
            else:
                pass
        pass
        

newLambda = np.array([])
HasLoop = 0
for i in range(len(n)):
    if i == 0:
        lamCalc =  2*n[0]*d[0] / m[0]
        newLambda = np.append(newLambda,lamCalc)
    elif i == (len(n)-1):
        lamCalc =  2*n[i]*d[i-1] / m[i+2]
        newLambda = np.append(newLambda,lamCalc)
    else:
        for k  in range(0,2):
            lamCalc =  2*n[i]*d[i-1+k] / m[i+k+HasLoop]
            print(n[i])
            print(d[i-1+k])
            print(m[i+k+HasLoop])
            print('')
            newLambda = np.append(newLambda,lamCalc)
            if k == 1: 
                HasLoop +=1
            else:
                pass