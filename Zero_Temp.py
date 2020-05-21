import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


#size of 1d lattice
lattice_size=10

#array of spins
spin_array=[1/2]

for i in np.arange(0,lattice_size,1):
    if spin_array[i] == 1/2 :
        spin_array.append((-1/2))
    else:
        spin_array.append((1/2))

print(spin_array)
#def hamiltonian(J,S1,S2):
    #random J
    #J = np.random.uniform(0,10)
    
