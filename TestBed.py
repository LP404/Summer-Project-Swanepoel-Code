import numpy as np
import Functions as F

#Factor in substrate
#Add in reflection?

#Mins

    # for k in range(len(vars()['nForMin'+str(i)])):
    #     if k == 0:
    #         dCalc =  (2*vars()['nForMin'+str(i)][0]) / (vars()['mForMin'+str(i)][0] * vars()['NewXForMin'+str(i)][0])
    #         vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
    #     elif k == (len(vars()['nForMin'+str(i)])-1):
    #         dCalc =  (2*vars()['nForMin'+str(i)][k]) / (vars()['mForMin'+str(i)][k] * vars()['NewXForMin'+str(i)][k])
    #         vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
    #     else:
    #         for l in range(0,2):
    #             lamCalc =  2*vars()['nForMin'+str(i)][k]* vars()['dForMin'+str(i)][k-1+l] / vars()['mForMin'+str(i)][[k+l+HasLoopMin]]
    #             vars()['NewXForMin'+str(i)] = np.append( vars()['NewXForMin'+str(i)],np.around(lamCalc,0))
            
    #         if l == 1: 
    #             HasLoopMin +=1
    #         else:
    #             pass  


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
            newLambda = np.append(newLambda,lamCalc)
            if k == 1: 
                HasLoop +=1
            else:
                pass
            
new_d = np.array([])
HasLoop1 = 0
for i in range(len(n)):
    if i == 0:
        dCalc =  (m[0] * newLambda[0]) / (2*n[0])
        new_d = np.append(new_d,dCalc)
    elif i == (len(n)-1):
        dCalc =  (m[i+2] * newLambda[i+2]) / (2*n[i])
        new_d = np.append(new_d,dCalc)
    else:
        for k  in range(0,2):
            dCalc =  (m[i+k+HasLoop1] * newLambda[i+k+HasLoop1]) / (2*n[k])
            print(i)
            print(i-1+k)
            print(i+k+HasLoop1)
            print('')
            new_d = np.append(new_d,dCalc)
            if k == 1: 
                HasLoop1 +=1
            else:
                pass
    
new_n = np.array([])

for i in range(len(newLambda)):
    nCalc = m[i]*newLambda[i] / (2*new_d[i])
    new_n = np.append(new_n,nCalc)

    
    
    for k in range(len(vars()['nForMax'+str(i)])):
        nCalc = vars()['mForMax'+str(i)][k]*vars()['NewXForMax'+str(i)][k] / (2*vars()['New_dForMax'+str(i)][k])
        vars()['New_nForMax'+str(i)] = np.append(vars()['New_nForMax'+str(i)],nCalc)

    for j in range(len(vars()['nForMin'+str(i)])):
        nCalc = vars()['mForMin'+str(i)][j]*vars()['NewXForMin'+str(i)][j] / (2*vars()['New_dForMin'+str(i)][j])
        vars()['New_nForMin'+str(i)] = np.append(vars()['New_nForMin'+str(i)],nCalc)