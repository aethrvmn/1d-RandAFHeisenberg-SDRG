import numpy as np

class Random_chain:

    version="v2.5"

    def __init__(self, number_of_bonds, ceiling, bottom):
        self.length = int(number_of_bonds)
        self.ceiling = ceiling
        self.bottom = bottom
        self.matrix_creator()
        self.bonds()
        self.parameters()
        self.rg_end()

    def matrix_creator(self): #creates and defines the matrix, adding "padding" of 2 zero rows+columns to fix the issue of biggest bond in the boundaries.
        self.bond_matrix = np.zeros(shape=(self.length ,self.length))
        initial_bonds = np.random.uniform(self.bottom, self.ceiling, self.length)
        np.fill_diagonal(self.bond_matrix, initial_bonds)
        A = np.c_[np.zeros(self.length), self.bond_matrix, np.zeros(self.length)]
        A = np.r_[[np.zeros(self.length + 2)], A, [np.zeros(self.length + 2)]] 
        self.bond_matrix = A
        return self

    def bonds(self): #finds the strongest bond and those next to it

        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()

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
        return self

    def parameters(self): #different physical parameters that are calculated
        self.distance = int(np.abs(self.right_index[1]-self.left_index[0]))
        self.system_energy = -(1/4)*np.sum(self.bond_matrix[self.bond_matrix > 0])
        self.gamma = -np.log(self.max_bond)
        return self

    def rescale(self): #not sure if I need this so I keep it here. If used, it should replace self.bonds in rg func.
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

    def renormalization(self): #func name says it all

        effective_bond = (self.left_bond*self.right_bond)/self.max_bond
        self.bond_matrix[self.left_index[0]][self.left_index[1]] = 0
        self.bond_matrix[self.right_index[0]][self.right_index[1]] = 0
        self.bond_matrix[self.left_index[0]][self.right_index[1]] = effective_bond
        self.bond_matrix[self.max_index[0]][self.max_index[1]] = self.max_bond - self.ceiling
        if self.bond_matrix[0][0]!=0:
            self.bond_matrix[0][0] = 0

        self.bonds()
        self.rg_end()
        return self

    def rg_end(self): #checks for non-singlets

        self.end_rg = 0
        self.bonds()
        if self.max_bond == 0.0:
            self.end_rg = 1
        return self
