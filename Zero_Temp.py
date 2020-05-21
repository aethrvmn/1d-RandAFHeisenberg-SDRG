import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


#size of 1d lattice
lattice_size=100

#array of spins
spin_array=[0.5]

#array of bond strenght
bonds_array=[]

#Creates the AF state
for i in np.arange(0,lattice_size,1):
    if spin_array[i] == 0.5 :
        spin_array.append((-0.5))
    else:
        spin_array.append((0.5))

for i in np.arange(0,lattice_size,1):
    J = np.random.uniform(0,10)
    bonds_array.append(J)



#finds the greatest bond
n=0
for i in np.arange(0,len(bonds_array),1):
    m = bonds_array[i]
    if m > n:
        x = bonds_array[i-1]
        y = bonds_array[i]
        z = bonds_array[i+1]
        n = y

greatest_bond=[x,y,z]

print(greatest_bond)

