import numpy as np

"Ok, let's try again, but this time let's make it work using matrices."

class ZT_Random_Spin:

    version="v0.1"

    def __init__(self, number_of_bonds, celing, floor):
        length = int(number_of_bonds) #the length of the chain (actually the total amount of bonds so chain-1)
        bond_matrix = numpy.zeros(shape=(length+2,length+2)) #initialzes the bond matrix
        initial_bonds = ceiling*np.random.rand(len(length))
        np.fill_diagonal(bond_matrix, initial_bonds) #adds the bonds
        self.average_strength() # compute the average strength of the bonds
        self.sys_energy() # compute the energy
        self.strongest_bond() #finds the strongest bond

    def sys_energy(self):
        self.system_energy = (-1/4)*np.sum(bond_matrix)
        return self

    def strongest_bond(self):
        self.max_bond = np.amax(bond_matrix)
        self.max_index = np.argmax(bond_matrix)
        return self

    def average_strength(self):
        self.mean = np.sum(bond_matrix)/length)
        return self

    def renormalization(self):
