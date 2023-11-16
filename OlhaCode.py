import numpy as np 
import cmath
import math
import matplotlib.pyplot as plt

pi = math.pi
Lam_array = np.array([188.280, 216.375, 237.397, 260.462, 299.328, 328.409, 377.414, 433.731, 475.871, 546.880, 628.485, 792.439, 953.900, 999.165])
epsilon_array = np.array([-2.64935707678035-0.607394806774214J, -3.48456652610259-0.916185620113791J, -4.58307563795461-1.26132338102961J, -6.0278895942636-1.58489319246111J, -7.57423617649664-2.0845305505626J, -9.96201307582504-2.74167851619151J, -12.5175882287977-3.77450064683337J, -16.4637561738609-4.96441144984948J, -20.687236496791-6.52944146772712J, -24.8336938204315-8.98915788187732J, -32.6625123372642-12.3754780289722J, -35.7864821001284-25.6990669252847J, -47.0681655960924-21.4081191092413J, -61.9063982422914-14.8559635938151J])

m = len(Lam_array)         
E = 30 #keV
v = math.sqrt((2*E*1000*1.6e-19)/(9.11e-31)) #m/s
beta = v/(3e+8)  
b = 361
theta_in = 0.123041
theta_out = 0.283794
theta_tilt = pi/4
#theta_in = pi*5/180    #for Polman's objective
#theta_out = pi/2
X_tilt = b/2*(1+2*theta_tilt/pi)     #coordinates of centre of reflecting objective (in terms of index)
Y_tilt = b/2


def I(m,theta,beta):

    lam = Lam_array[m]
    eps = epsilon_array[m]
    K = (((1/137)*beta**2) / (pi*lam))          #coefficients 
    S = math.sin(theta)**2 *math.cos(theta)**2
    Anum = abs(((eps - 1)*(1 - beta**2 + beta*cmath.sqrt(eps - (math.sin(theta)**2)))))**2
    Aden = abs(((1 - beta**2*math.cos(theta)**2)* (1 + beta*cmath.sqrt((eps - (math.sin(theta))**2))) * ((eps*math.cos(theta)) + cmath.sqrt(eps - (math.sin(theta))**2))))**2  
    A = Anum/Aden
    Ires = K*S*A
    return Ires
    
I_array = np.zeros((b,b))
Theta_array = np.zeros((b,b))
Phi_array = np.zeros((b,b))
Omega_array = np.zeros((b,b))
sin = math.sin(pi/(2*b))
for i in range(b):
    for j in range(b):

        phi = np.arctan2((j-b/2),(i-b/2))
        theta = pi/b*math.sqrt((i-b/2+1/2)**2 + (j-b/2+1/2)**2) 
        Theta_array[i,j] = theta
        Phi_array[i,j] = phi
        if theta <= pi/2:
            Imean = 0
            for k in range(m):
                Imean += I(k,theta, beta)
            I_array[i,j] = (Imean/m)  
            if theta > 0:
                Omega_array[i,j]= 2*pi*math.sin(Theta_array[i,j])*sin/(Theta_array[i,j]*b)  #solid angle per pixel
            else:
                Omega_array[i,j]= 2*pi*sin/b

Sum = np.sum(Omega_array)
print(Sum)


Mask_array = np.zeros((b,b))
for i in range(b):
    for j in range(b):
        x = math.sin(Theta_array[i,j])*math.cos(Phi_array[i,j])       #spherical->cartesian
        y = math.sin(Theta_array[i,j])*math.sin(Phi_array[i,j])
        z = math.cos(Theta_array[i,j])
        
        x2 = x                                                      #rotation Rx
        y2 = math.cos(theta_tilt)*y-math.sin(theta_tilt)*z
        z2 = math.sin(theta_tilt)*y+math.cos(theta_tilt)*z
       
        theta_test = np.arctan2(math.sqrt(x2**2+y2**2),z2)            #back to the cartesian, theta only
        
        if theta_test < theta_out:
            if theta_test > theta_in:
                Mask_array[i,j] = 1 
        
        #Mask_array[i,j] = theta_test
        
C = np.multiply(I_array, Mask_array)       #element-wise multiplication , collected intensity 
C = np.multiply(C, Omega_array)
s = np.sum(C)                                           #sum of the collected intensity per nm
print(s)


fig, ax1 = plt.subplots()
im1 = ax1.imshow(Theta_array, cmap='viridis')    #or plasma
plt.colorbar(im1, ax = ax1)
fig, ax2 = plt.subplots()
im2 = ax2.imshow(I_array, cmap='viridis')   #in angle space, hemisphere, 2D map for fixed wavelength
plt.colorbar(im2, ax = ax2)
fig, ax3 = plt.subplots()
im3 = ax3.imshow(Mask_array, cmap='viridis')  
ax3.set_title('Mask', fontsize = 20)
plt.colorbar(im3, ax = ax3)
fig, ax4 = plt.subplots()
im4 = ax4.imshow(C, cmap='viridis')
ax4.set_title('Mask*intensity map', fontsize = 20)
plt.colorbar(im4, ax = ax4)
fig, ax5 = plt.subplots()
im5 = ax5.imshow(Phi_array, cmap='viridis')    #or plasma
plt.colorbar(im5, ax = ax5)
fig, ax6 = plt.subplots()
im6 = ax6.imshow(Omega_array, cmap='viridis')    #or plasma
plt.colorbar(im6, ax = ax6)
