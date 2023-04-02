import numpy as np

# np.random.seed(20) # Comment if you want randomised results, it was used to get the graphs for the thesis text.


class Random_chain:

    def __init__(self, bond_number, resc = True):
        self.length = int(bond_number)
        self.ceiling = 1
        self.matrix_creator()
        self.bonds()
        self.resc = resc

    def matrix_creator(self):
        initial_bonds = np.random.uniform(0, self.ceiling, self.length)
        padded_length = self.length + 2
        self.bond_matrix = np.zeros((padded_length, padded_length))
        np.fill_diagonal(self.bond_matrix[1:-1, 1:-1], initial_bonds)
        self.val = np.amax(self.bond_matrix)
        return self

    def bonds(self):
        self.max_bond = np.amax(self.bond_matrix)
        self.max_index = np.argwhere(self.bond_matrix.max() == self.bond_matrix).ravel()

        if self.max_bond != 0:
            search_x, search_y = self.max_index

            left_bonds = self.bond_matrix[:search_x, :search_y]
            left_indices = np.argwhere(left_bonds > 0)
            if left_indices.size > 0:
                self.left_index = left_indices[-1]
            else:
                self.left_index = np.array([0, 0])
            self.left_bond = self.bond_matrix[tuple(self.left_index)]

            right_bonds = self.bond_matrix[search_x + 1:, search_y + 1:]
            right_indices = np.argwhere(right_bonds > 0)
            if right_indices.size > 0:
                self.right_index = right_indices[0] + [search_x + 1, search_y + 1]
            else:
                self.right_index = np.array([self.length + 1, self.length + 1])
            self.right_bond = self.bond_matrix[tuple(self.right_index)]

            self.eff_bond = self.right_bond * self.left_bond / self.max_bond
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
        return self

    def renormalization(self): # Self-explanatory
        self.decimate()
        if self.resc:
            self.rescale()
        self.bonds()
        return self
