import numpy as np

Array = np.array([[1,1,1,4,2,2,3,3,3],[np.pi,10,4,7,-1,4,7,18,9]])

# def DYThorLabs(Array):
    
newXarray = np.array([])
newYarray = np.array([])

newXarray = np.append(newXarray,Array[0][0])
backTotalCount = 0
backTotal = 0

for i in range(1,len(Array[0])):
    if Array[0][i] == Array[0][i-1]:
        backTotalCount += 1
        
    else:   
        newXarray = np.append(newXarray,Array[0][i])
        newYarray = np.append(newYarray,np.mean(Array[1][backTotal:i]))
        backTotal = backTotalCount + 1


backTotal = backTotalCount + 1
newYarray = np.append(newYarray,np.mean(Array[1][backTotal:i]))

FixedArray = np.array([newXarray,newYarray])
    
    # return FixedArray


