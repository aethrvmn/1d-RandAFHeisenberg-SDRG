import numpy as np

"Ok, let's try again, but this time let's make it work using matrices."

class ZT_Random_Spin:

    version="v0.9.7"

    def __init__(self, number_of_bonds, ceiling, floor):
        self.length = int(number_of_bonds) #the self.length of the chain (actually the total amount of bonds so chain-1)
        self.ceiling = ceiling
        self.floor = floor
        self.matrix_creator()
        self.average_strength()
        self.bonds()
        self.sys_energy()
        self.logarithmic()


    #creates and defines the matrix, adding "padding" of 2 zero rows+columns to dodge the issue of biggest bond in the boundaries
    def matrix_creator(self):
        self.bond_matrix = np.zeros(shape=(self.length ,self.length)) #initialzes the bond matrix
        initial_bonds = self.ceiling*np.random.uniform(0, self.ceiling, self.length)
        np.fill_diagonal(self.bond_matrix, initial_bonds) #adds the bonds
        A = np.c_[np.zeros(self.length), self.bond_matrix, np.zeros(self.length)] #this allows us to add two columns of 0s, essentially closing the system and also solves IndexErrors for when the max bond is in the boundaries
        A = np.r_[[np.zeros(self.length + 2)], A, [np.zeros(self.length + 2)]] #this adds the rows of 0s, see above
        A = np.c_[np.zeros(self.length + 2), A, np.zeros(self.length + 2)]
        A = np.r_[[np.zeros(self.length + 4)], A, [np.zeros(self.length + 4)]]
        A = np.c_[np.zeros(self.length + 4), A, np.zeros(self.length + 4)]
        A = np.r_[[np.zeros(self.length + 6)], A, [np.zeros(self.length + 6)]]
        self.bond_matrix = A

    #computes the energy
    def sys_energy(self):
        self.system_energy = (-1/4)*np.sum(self.bond_matrix)
        return self

    #finds the strongest bond and those next to it
    def bonds(self):
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()
        self.left_index = np.array([self.max_index[0] - 1, int(np.nonzero(self.bond_matrix[][self.max_index[0] - 1]))])
        print(self.left_index)
        self.right_index = np.array([self.max_index[1] + 1, int(np.nonzero(self.bond_matrix[self.max_index[1] + 1])[0])])
        self.left_bond = self.bond_matrix[self.left_index[0]][self.left_index[1]]
        self.right_bond = self.bond_matrix[self.right_index[0]][self.right_index[1]]
        try:
        print(self.left_index, self.max_index, self.right_index)
        except TypeError:
        print('daisy')
        return self

    #computes the average strength of the bonds
    def average_strength(self):
        self.mean = np.sum(self.bond_matrix)/self.length
        return self

    def logarithmic(self):
        self.logmax = -np.log(self.max_bond)

    # def RG_logic(self):
    #     #self.bond_matrix[self.bond_matrix == 0.0] = np.nan
    #     self.new_rg_spot = self.bond_matrix[self.left_index[0]][self.left_index[0]] != 0 and self.bond_matrix[self.right_index[0]][self.right_index[0]] != 0
    #     self.right_rg_spot = self.bond_matrix[self.left_index][self.left_index] != 0 and self.bond_matrix[self.right_index][self.right_index] == 0
    #     self.left_rg_spot = self.bond_matrix[self.left_index][self.left_index] == 0 and self.bond_matrix[self.right_index][self.right_index] != 0
    #     self.old_rg_spot = self.max_index[0] != self.max_index[1]# and
    #     return self

    #This is the RG process
    #This thing is extremely hard, I have no hopes but here it goes
    def renormalization(self):
        while (self.max_bond > self.floor):
            # self.RG_logic()
            energy_prime = -(1/4)*self.max_bond - ((3/16)/self.max_bond)*((self.left_bond**2) + (self.right_bond**2)) # find the energy contribution
            bond_prime = (self.left_bond*self.right_bond)/(2*self.max_bond) # find the strength of the new bond that will exist after we remove the spins
            self.local_energy = (-1/4)*np.sum(self.bond_matrix[self.max_index])
            self.new_local_energy = energy_prime -(1/4)*bond_prime # finds the energy that we will add to the total energy after we remove the spins/bonds that existed
            self.bond_matrix[self.right_index[0]] = 0.00000000001
            self.bond_matrix[self.left_index[0]] = 0.00000000001
            self.bond_matrix[self.left_index[0]][self.right_index[1]] = bond_prime
            #print(self.left_index, self.right_index)
        #This is the first order of RG
            #This is for the situation where the next bond has been RG transformed.
            # elif self.right_rg_spot:
            #     self.bond_matrix[self.left_index[0]][self.right_index[1]] = bond_prime
            #     print(self.left_index, self.right_index, 'nrg2')
            # #This is for the situation where the previous bond has been RG transformed.
            # elif self.left_rg_spot:
            #     self.bond_matrix[self.left_index[0]][self.right_index[1]] = bond_prime
            #     print(self.left_index, self.right_index,'nrg3')
            #This is for the situation where a bond that has already been transformed is being RG transformed again.
            #elif self.old_rg_spot:
                # self.bond_matrix[self.left_index][self.right_index + 1] = bond_prime
                # self.bond_matrix[self.right_index + 1] = 0
            #print('WIP')
            self.bond_matrix[self.max_index[0]][self.max_index[1]] = self.max_bond - self.ceiling
            self.system_energy = self.system_energy - self.local_energy + self.new_local_energy # calculates the new energy of the chain by removing the previous contribution of the strongest bond and adding the new contribution of the newly weak bond in the same spot
            self.average_strength() # recalculates the average strength
            self.logarithmic()
            self.bonds() # refreshes the strongest bond
        return self
