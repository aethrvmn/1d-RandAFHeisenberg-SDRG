import numpy as np


'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain and the RG transformation (elimination transformation) for a non-zero temperature system. Currently there is an issue with what happens when the size of the chain falls to zero, and shows an error. I also need to change the distribution method for the bonds so it's an actual random distribution to properly study the RG flow.
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

    # Calculates the mean value of the bonds
    def mean_bonds(self):
        self.mean = np.mean(self.bonds)
    
    # This finds the strongest bond
    def strong_bond(self):
        self.mega_bond = np.amax(self.bonds) # strongest bond
        self.mega_index = np.argmax(self.bonds) # the position of the strongest bond
        self.mega_state = self.state[self.mega_index] # the state of the strongest bond

        try:
            self.left_state = self.state[self.mega_index - 1]
        except IndexError:
            self.left_state = 0
        try:
            self.right_state = self.state[self.mega_index + 1]
        except IndexError:
            self.right_state = 0

        try:
            self.left_mini_bond = self.bonds[self.mega_index-1] # the bond between the left spin and the one left to it
        except IndexError:
            self.left_mini_bond = 0
        try:
            self.right_mini_bond = self.bonds[self.mega_index+1] # the bond between the right spin and the one to the right of it
        except IndexError:
            self.right_mini_bond = 0
        self.local_hamiltonian = self.mega_bond*self.state[self.mega_index]*self.right_state # finds the energy that we will remove from the total energy during the transformation
        return self
    
    #Equation of the free energy of the non zero temprature system
    def sys_free_energy(self):
        self.free_zero = -0.75*self.mega_bond - (1/self.beta)*np.log(1+(3*np.exp(-self.beta*self.mega_bond)))
        
        return self
    
    
    
    
    
    
    #We calculate the Hamiltonian, which in the zero temp case would be the system energy
    #def hamiltonian(self):
        #spin_states1 = np.delete(self.state, [0]) # creates one of the two pseudo-chains to help compute the energy
        #spin_states2 = np.delete(self.state, [self.length-1]) # this is the second pseudo-chain
        #self.hamiltonian_energy = np.sum(self.bonds*spin_states1*spin_states2)
        #return self
    
    #We calculate the partition function to help with the evaluation of the free energy
    #def partition_function(self):
        #self.hamiltonian()
        #self.exponential = np.exp(-self.beta*self.hamiltonian_energy)
        #self.partition = np.sum(self.exponential)
        #return self

    #Instead of the ground energy, since we have a finite temperature, we calculate the free energy
    #def free_energy(self):
        #self.partition_function()
        #self.energy = - (1/self.beta)*np.log(self.partition)
        #return self
