import numpy as np

"Ok, let's try again, but this time let's make it work using matrices."

class ZT_Random_Spin:

    version="v0.5"

    def __init__(self, number_of_bonds, ceiling, floor):
        self.length = int(number_of_bonds) #the self.length of the chain (actually the total amount of bonds so chain-1)
        self.bond_matrix = np.zeros(shape=(self.length+2,self.length+2)) #initialzes the bond matrix
        initial_bonds = ceiling*np.random.rand(self.length)
        np.fill_diagonal(self.bond_matrix, initial_bonds) #adds the bonds
        self.average_strength() # compute the average strength of the bonds
        self.sys_energy() # compute the energy
        self.strongest_bond() #finds the strongest bond

    def sys_energy(self):
        self.system_energy = (-1/4)*np.sum(self.bond_matrix)
        return self

    def strongest_bond(self):
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix)
        return self

    def average_strength(self):
        self.mean = np.sum(self.bond_matrix)/self.length
        return self

    #def renormalization(self):
