import numpy as np
import scipy.special as sp
"Ok, let's try again, but this time let's make it work using matrices."

class NZT_Random_Spin:

    version="v0.1"

    def __init__(self, number_of_bonds, ceiling, floor, temperature):
        self.length = int(number_of_bonds) #the self.length of the chain (actually the total amount of bonds so chain-1)
        self.bond_matrix = np.zeros(shape=(self.length+2,self.length+2)) #initialzes the bond matrix
        initial_bonds = ceiling*np.random.rand(self.length)
        self.beta = 1/temperature
        np.fill_diagonal(self.bond_matrix, initial_bonds) #adds the bonds
        self.average_strength() # compute the average strength of the bonds
        self.sys_free_energy() # compute the energy
        self.strongest_bond() #finds the strongest bond
        self.chi = self.beta*self.max_bond


    def sys_free_energy(self):
        hamiltonian = (-1/4)*self.bond_matrix
        self.free_energy = (-1/self.beta)*sp.logsumexp(-self.beta*hamiltonian)
        return self

    def strongest_bond(self):
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix)
        return self

    def average_strength(self):
        self.mean = np.sum(self.bond_matrix)/self.length
        return self

    def factor_functions(self):
        self.V = (1-(np.exp(-self.chi))*(1-self.chi))/(1+(3*np.exp(-self.chi)))
        self.W = (1-(np.exp(-self.chi))*(1+self.chi))/(1+(3*np.exp(-self.chi)))
