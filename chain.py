import numpy as np

# I'm sure this joke will age well
np.random.seed(69420) # Comment if you want randomised results, it was used to get the graphs for the thesis text.
# I'm sure this joke will age well

class Random_chain:

    version = "v5"

    def __init__(self, bond_number, ceiling):
        self.length = int(bond_number)
        self.ceiling = ceiling
        self.matrix_creator()
        self.bonds()

    def matrix_creator(self): # Creates and defines the matrix, adding "padding" of 2 zero rows+columns to fix the issue of biggest bond in the boundaries.
        matrix = np.zeros(shape=(self.length ,self.length))
        initial_bonds = np.random.uniform(0, self.ceiling, self.length)
        np.fill_diagonal(matrix, initial_bonds)
        matrix = np.c_[np.zeros(self.length), matrix, np.zeros(self.length)]
        matrix = np.r_[[np.zeros(self.length + 2)], matrix, [np.zeros(self.length + 2)]] 
        self.bond_matrix = matrix
        self.val = np.amax(self.bond_matrix)
        return self

    def bonds(self): # Searches for the strongest bond, those next to it, and also converts them to their log counterparts
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()
        
        if self.max_bond != 0:
            search_x = self.max_index[0]
            search_y = self.max_index[1]

            self.left_index = np.array([0,0])
            for i in range(search_x):
                for j in range(search_y):
                    if self.bond_matrix[i][j] > 0:
                        self.left_index = np.array([i, j])
            self.left_bond = self.bond_matrix[self.left_index[0]][self.left_index[1]]
            
            self.right_index = np.array([self.length+1,self.length+1])
            for i in range(self.length, search_x, -1):
                for j in range(self.length, search_y, -1):
                    if self.bond_matrix[i][j] > 0:
                        self.right_index = np.array([i, j])
            self.right_bond = self.bond_matrix[self.right_index[0]][self.right_index[1]]

            self.eff_bond = self.right_bond*self.left_bond/(self.max_bond)

            if self.eff_bond == 0:
                self.zeta = 0
            else:
                self.zeta = np.log(self.max_bond/self.eff_bond)
            self.gamma = -np.log(self.max_bond)
            self.eta = self.zeta/self.gamma
        return self

    def decimate(self): # First Step of the RG
        self.bond_matrix[self.max_index[0]][self.max_index[1]] = -self.max_bond
        self.bond_matrix[self.left_index[0]][self.left_index[1]] = 0
        self.bond_matrix[self.right_index[0]][self.right_index[1]] = 0
        self.bond_matrix[self.left_index[0]][self.right_index[1]] = self.eff_bond
        return self

    def rescale(self): # Second Step of the RG
        self.max_bond = np.amax(self.bond_matrix)
        if self.max_bond != self.val:
            self.bond_matrix[self.bond_matrix > 0] += np.abs(self.val-self.max_bond)
        self.bonds()
        return self

    def renormalization(self): # Self-explanatory
        self.decimate()
        self.rescale()
        return self
    
    def parameters(self): # Different physical parameters that can be calculated
        self.distance = int(np.abs(self.right_index[1]-self.left_index[0]))
        self.system_energy = -(1/4)*np.sum(self.bond_matrix[self.bond_matrix > 0])
        return self