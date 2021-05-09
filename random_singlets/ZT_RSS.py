import numpy as np
#%%
class ZT_Random_Spin:

    version="v2.1"

    def __init__(self, number_of_bonds, ceiling):
        self.length = int(number_of_bonds) #the length of the chain (actually the total amount of bonds so chain-1)
        self.ceiling = ceiling
        self.matrix_creator()
        self.average_strength()
        self.bonds()
        self.rg_end()
        self.sys_energy()
        self.logarithmic()

    #creates and defines the matrix, adding "padding" of 2 zero rows+columns to dodge the issue of biggest bond in the boundaries
    def matrix_creator(self):
        self.bond_matrix = np.zeros(shape=(self.length ,self.length)) #initialzes the bond matrix
        initial_bonds = self.ceiling*np.random.uniform(0, self.ceiling, self.length)
        np.fill_diagonal(self.bond_matrix, initial_bonds) #adds the bonds
        A = np.c_[np.zeros(self.length), self.bond_matrix, np.zeros(self.length)] #this allows us to add two columns of 0s, essentially closing the system and also solves IndexErrors for when the max bond is in the boundaries
        A = np.r_[[np.zeros(self.length + 2)], A, [np.zeros(self.length + 2)]] #this adds the rows of 0s, see above
        self.bond_matrix = A
        return self

    #computes the energy
    def sys_energy(self):
        self.system_energy = -(1/4)*np.sum(self.bond_matrix[self.bond_matrix > 0])
        return self

    #finds the strongest bond and those next to it
    def bonds(self):
        #finds the strongest bond
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()

        #finds the bond to the left
        self.left_index = np.array([0, 0])
        self.left_bond = 0
        left_search_area_x = self.max_index[0]
        left_search_area_y = self.max_index[1]
        if self.max_index[0] != 1:
            for i in np.arange(left_search_area_x):
                for j in np.arange(left_search_area_y):
                    if self.bond_matrix[i][j] > 0:
                        self.left_index = np.array([i, j])
                        self.left_bond = self.bond_matrix[i][j]
        else:
            self.left_index = np.array([0, 0])
            self.left_bond = 0

        #finds the bond to the right
        self.right_index = np.array([self.length+1, self.length+1])
        self.right_bond =0
        if self.max_index[0] != self.length:
            for i in np.arange(self.length, self.max_index[0], -1):
                for j in np.arange(self.length, self.max_index[1], -1):
                    if self.bond_matrix[i][j] > 0:
                        self.right_index = np.array([i, j])
                        self.right_bond = self.bond_matrix[i][j]
        else:
            self.right_index = np.array([self.length+1, self.length+1])
            self.right_bond =0
        #print(self.left_index, self.max_index, self.right_index)
        return self

    #computes the average strength of the bonds
    def average_strength(self):
        self.mean = np.sum(self.bond_matrix)/self.length
        return self

    def logarithmic(self):
        self.logmax = -np.log(self.max_bond)
        return self

    def rgflow(self):
        self.flow = self.bond_matrix[np.nonzero(self.bond_matrix)]
        return self
    
    
    def rescale(self):
        temporary_val = self.max_bond
        self.bonds()
        val = temporary_val - self.max_bond
        for i in np.arange(self.length+2):
            for j in np.arange(self.length+2):
                if self.bond_matrix[i][j] > 0:
                    self.bond_matrix[i][j] += val
                else:
                    continue
        return self

    #This is the RG process
    def renormalization(self):
        self.bond_matrix[self.left_index[0]][self.left_index[1]] = 0
        self.bond_matrix[self.right_index[0]][self.right_index[1]] = 0
        self.bond_matrix[self.left_index[0]][self.right_index[1]] = bond_prime
        self.bond_matrix[self.max_index[0]][self.max_index[1]] = self.max_bond - self.ceiling
        if self.bond_matrix[0][0]!=0:
            self.bond_matrix[0][0] = 0

        self.distance = self.right_index[1]-self.left_index[0]
        self.average_strength() # recalculates the average strength
        self.sys_energy()
        #self.logarithmic() #evaluates the log variable for the rg flow
        self.rescale() # rescales the matix and also refreshes the biggest bond search
        self.rg_end() #checks for non-singlets
        return self

    def rg_end(self):
        self.end_rg = 0
        self.bonds()
        if self.max_bond == 0.0:
            self.end_rg = 1
        return self

# %%
