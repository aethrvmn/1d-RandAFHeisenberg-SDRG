import numpy as np

"Ok, let's try again, but this time let's make it work using matrices."

class ZT_Random_Spin:

    version="v0.1"

    def __init__(self, number_of_bonds, celing, floor):
        length = int(number_of_bonds) #the length of the chain (actually the total amount of bonds so chain-1)
        self.bond_matrix = numpy.zeros(shape=(self.length,self.length)) #initialzes the bond matrix
        self.max_bond = 0 #initialzes a value for the biggest bond
        #this creates the bond matrix
        for i in range(length):
            for j in range(length):
                if i == j:
                    self.bond_matrix[i][j] = ceiling*np.random.uniform(0,1)
                    if self.max_bond < self.bond_matrix[i][j]:
                        self.max_bond = self.bond_matrix[i][j]
                    else:
                        continue
                else:
                    continue
        #end of bond_matrix
        self.average_strength() # compute the mean
        self.sys_energy() # compute the energy

    def sys_energy(self):
        energy = np.zeros(shape=(int(length)))
        for i in range(length):
            for j in range(length):
                energy[i] += (-1/4)*bond_matrix[i][j]
        system_energy = np.sum(energy)

    def average_strength(self):
        mean = np.zeros(shape=(int(length)))
        for i in range(length):
            for j in range(length):
                mean[i] += bond_matrix[i][j]
         average = np.sum((mean/length))
