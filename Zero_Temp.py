import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from classes.spin_chain import spin_chain

#Size of 1D Lattice
lattice_size=1000

#Array of Spins
spin_array=[0.5]

#Array of Bond Strenght
bonds_array=[]

#Creates the AF state
for i in np.arange(0,lattice_size,1):
    if spin_array[i] == 0.5 :
        spin_array.append((-0.5))
    else:
        spin_array.append((0.5))
#Creates Randomly
for i in np.arange(0,lattice_size,1):
    J = np.random.normal(0,5)
    if J < 0:
        J = -J
    bonds_array.append(J)

#Finds the Greatest Bond
n=0
for i in np.arange(0,len(bonds_array),1):
    m = bonds_array[i]
    if m > n:
        x = bonds_array[i-1]
        y = bonds_array[i]
        z = bonds_array[i+1]
        n = y
    plt.scatter(i,m)
    plt.pause(0.1)

greatest_bond=[x,y,z]
plt.show()
print(greatest_bond)
