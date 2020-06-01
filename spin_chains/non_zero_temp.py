import numpy as np
import scipy.special as sp

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
        self.beta = 1/temperature
        #self.floor = floor
        self.mean_bonds() # compute the mean
        self.strong_bond()
        self.chi = self.beta*self.mega_bond
        self.sys_free_energy()

    # Calculates the mean value of the bonds
    def mean_bonds(self):
        self.mean = np.mean(self.bonds)

    # This finds the strongest bond
    def strong_bond(self):
        self.mega_bond = np.amax(self.bonds) # strongest bond
        self.mega_index = np.argmax(self.bonds) # the position of the strongest bond
        try:
            self.mega_state = self.state[self.mega_index] # the state of the strongest bond
        except IndexError:
            print('No more bonds left')
        #This finds the states of the spins next to the ones sharing the bond
        try:
            self.left_state = self.state[self.mega_index - 1]
        except FutureWarning:
            self.left_state = 0
            
        try:
            self.second_spin = self.state[self.mega_index+1]
        except IndexError:
            self.second_spin = 0
            
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
        
        self.local_hamiltonian = self.mega_bond*self.state[self.mega_index]*self.second_spin # finds the energy that we will remove from the total energy during the transformation
        self.local_free_energy = -(1/self.beta)*np.log(np.exp(-self.beta*self.local_hamiltonian)) # in essence the local hamiltonian
        return self

    #Here we define the V and W functions that are used to compute the free energy
    def factor_functions(self):
        self.V = (1-(np.exp(-self.chi))*(1-self.chi))/(1+(3*np.exp(-self.chi)))
        self.W = (1-(np.exp(-self.chi))*(1+self.chi))/(1+(3*np.exp(-self.chi)))

    #Equation of the free energy of the non zero temprature system
    def sys_free_energy(self):
        spin_states1 = np.delete(self.state, [0]) # creates one of the two pseudo-chains to help compute the energy
        spin_states2 = np.delete(self.state, [self.length-1]) # this is the second pseudo-chain
        self.hamiltonian = np.sum(self.bonds*spin_states1*spin_states2)
        self.free_zero = -0.75*self.mega_bond - (1/self.beta)*np.log(1+(3*np.exp(-self.beta*self.mega_bond)))
        self.free_energy = (-1/self.beta)*sp.logsumexp(-self.beta*self.hamiltonian)
        return self
    # This is the RG process
    
    def nzt_elimination(self):
        # while (self.mega_bond > self.floor):
        self.factor_functions()
        self.free_prime = self.free_zero - (0.1875*(self.left_mini_bond**2 + self.right_mini_bond**2)*self.V)
        try:
            self.bond_prime = (self.left_mini_bond*self.right_mini_bond*self.W)/2*self.mega_bond
        except IndexError:
            self.bond_prime = 0
        self.exponential = np.sum(np.exp(-self.beta*self.bond_prime*self.left_state*self.right_state))
        self.new_local_free_energy = self.free_prime - (1/self.beta)*np.log(self.exponential)
        
        self.state = np.delete(self.state, [self.mega_index, self.mega_index+1]) # the new chain after removing the spins sharing the strongest bond
        self.bonds[self.mega_index]=self.bond_prime # replacing the strongest bond with the bond prime
        self.bonds = np.delete(self.bonds, [self.mega_index-1, self.mega_index+1]) # removing the bonds next to the bond prime to get the proper new chain
        
        self.free_energy = self.free_energy - self.local_free_energy + self.new_local_free_energy
        self.strong_bond() # refreshes the strongest bond
        self.mean_bonds() # finds the new mean bond
        #return self
