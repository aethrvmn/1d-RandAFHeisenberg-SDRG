import numpy as np


'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain
and the RG transformation (elimination transformation) for a non-zero temperature system.
Currently there is an issue with what happens when the size of the chain falls to zero,
and shows an error.

I also need to change the distribution method for the bonds so it's an actual random distribution
to properly study the RG flow.
'''
class NZT_Chain:

    version = 'v0.1'


    # Initial call, starts the project by taking in the length and the maximum energy ceiling, and creating a random spin chain based on those.
    def __init__(self, length, ceiling, temperature):#, floor):
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.bonds = ceiling*np.random.rand(len(self.state) - 1)
        self.length = len(self.state)
        #self.floor = floor
        self.mean_bonds() # compute the mean
        self.strong_bond()
        self.beta = 1/temperature
        self.chi = self.beta*self.mega_bond
        self.sys_free_energy()

    # Calculates the mean value of the bonds
    def mean_bonds(self):
        self.mean = np.mean(self.bonds)

    # This finds the strongest bond
    def strong_bond(self):
        self.mega_bond = np.amax(self.bonds) # strongest bond
        self.mega_index = np.argmax(self.bonds) # the position of the strongest bond
        self.mega_state = self.state[self.mega_index] # the state of the strongest bond
        #This finds the states of the spins next to the ones sharing the bond
        try:
            self.left_state = self.state[self.mega_index - 1]
        except IndexError:
            self.left_state = 0
        try:
            self.right_state = self.state[self.mega_index + 2]
        except IndexError:
            self.right_state = 0

        try:
            self.left_mini_bond = self.bonds[self.mega_index-1] # the bond between the left spin and the one left to it
        except IndexError:
            self.left_mini_bond = 0
        try:
            self.right_mini_bond = self.bonds[self.mega_index+2] # the bond between the right spin and the one to the right of it
        except IndexError:
            self.right_mini_bond = 0
        self.local_hamiltonian = self.mega_bond*self.state[self.mega_index]*self.right_state # finds the energy that we will remove from the total energy during the transformation
        return self

    #Here we define the V and W functions that are used to compute the free energy
    def factor_functions(self):
        self.V = (1-(np.exp(-self.chi))*(1-self.chi))/(1+(3*np.exp(-self.chi)))
        self.W = (1-(np.exp(-self.chi))*(1+self.chi))/(1+(3*np.exp(-self.chi)))

    #Equation of the free energy of the non zero temprature system
    def sys_free_energy(self):
        self.factor_functions()
        self.free_zero = -0.75*self.mega_bond - (1/self.beta)*np.log(1+(3*np.exp(-self.beta*self.mega_bond)))
        self.free_prime = self.free_zero - (0.1875*(self.left_mini_bond**2 + self.right_mini_bond**2)*self.V)
        self.bond_prime = (self.left_mini_bond*self.right_mini_bond*self.W)/2*self.mega_bond
        self.exponential = np.sum(np.exp(-self.beta*self.bond_prime*self.left_state*self.right_state))
        self.free_energy = self.free_prime - (1/self.beta)*np.log(self.exponential)
        return self
    # This is the RG process
    def nzt_elimination_transformation(self):
       # while (self.mega_bond > self.floor):
       

        self.strong_bond() # refreshes the strongest bond
        self.mean_bonds()
        return self
