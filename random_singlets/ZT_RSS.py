import numpy as np

"Ok, let's try again, but this time let's make it work using matrices."

class ZT_Random_Spin:

    version="v0.9"

    def __init__(self, number_of_bonds, ceiling, floor):
        self.length = int(number_of_bonds) #the self.length of the chain (actually the total amount of bonds so chain-1)
        self.bond_matrix = np.zeros(shape=(self.length,self.length)) #initialzes the bond matrix
        initial_bonds = ceiling*np.random.rand(self.length)
        np.fill_diagonal(self.bond_matrix, initial_bonds) #adds the bonds
        A = np.c_[np.zeros(self.length), self.bond_matrix, np.zeros(self.length)] #this allows us to add two columns of 0s, essentially closing the system and also solves IndexErrors for when the max bond is in the boundaries
        self.bond_matrix = np.r_[[np.zeros(self.length+2)], A, [np.zeros(self.length+2)]] #this adds the rows of 0s, see above
        self.average_strength()
        self.bonds()
        self.sys_energy()
        self.logarithmic()
        self.floor = floor


    #computes the energy
    def sys_energy(self):
        self.system_energy = (-1/4)*np.sum(self.bond_matrix)
        return self

    #finds the strongest bond and those next to it
    def bonds(self):
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()
        self.left_index = self.max_index[0]-1
        self.right_index = self.max_index[0]+1
        self.left_bond = np.sum(self.bond_matrix[self.left_index])
        self.right_bond = np.sum(self.bond_matrix[self.right_index])
        return self


    #computes the average strength of the bonds
    def average_strength(self):
        self.mean = np.sum(self.bond_matrix)/self.length
        return self

    def logarithmic(self):
        self.logmax = -np.log(self.max_bond)
        self.logbonds = np.log(self.max_bond/self.bond_matrix)
        self.logbonds[self.logbonds == np.inf] = 0
    #This is the RG process
    def renormalization(self):
        #while (self.max_bond > self.floor):
        energy_prime = -(1/4)*self.max_bond - ((3/16)/self.max_bond)*((self.left_bond**2) + (self.right_bond**2)) # find the energy contribution
        bond_prime = (self.left_bond*self.right_bond)/(2*self.max_bond) # find the strength of the new bond that will exist after we remove the spins
        self.local_energy = (-1/4)*np.sum(self.bond_matrix[self.max_index])
        self.new_local_energy = energy_prime -(1/4)*bond_prime # finds the energy that we will add to the total energy after we remove the spins/bonds that existed

        # self.state = np.delete(self.state, [self.mega_index, self.mega_index+1]) # the new chain after removing the spins sharing the strongest bond
        # self.bonds[self.mega_index]=self.bond_prime # replacing the strongest bond with the bond prime
        # self.bonds = np.delete(self.bonds, [self.mega_index-1, self.mega_index+1]) # removing the bonds next to the bond prime to get the proper new chain

        self.system_energy = self.system_energy - self.local_energy + self.new_local_energy # calculates the new energy of the chain by removing the previous contribution of the strongest bond and adding the new contribution of the newly weak bond in the same spot
        self.average_strength() # recalculates the average strength
        self.bonds() # refreshes the strongest bond
        self.logarithmic()
        return self
