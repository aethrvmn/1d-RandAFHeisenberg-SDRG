import numpy as np
import matplotlib.pyplot as plt



'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain and the RG transformation (elimination transformation) for a zero temprature system. Currently there is an issue with what happens when the size of the chain falls to zero, and shows an error. I also need to change the distribution method for the bonds so it's an actual random distribution to properly study the RG flow.
'''
class Chain:
    
    version='v1.1'
    

    #Initial call, starts the project by taking in the length and the maximum energy ceiling, and creating a random spin chain based on those.
    def __init__(self, length, ceiling):
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.bonds = ceiling*np.random.rand(len(self.state) - 1)
        self.length = len(self.state)
        self.sys_energy() # compute the energy
        self.strong_bond() # find the strongest bond
        self.mean_bonds()
        
    #Calculates the mean value of the bonds    
    def mean_bonds(self):
        self.mean = np.mean(self.bonds)
    
    #This calculates the system energy    
    def sys_energy(self):
        
        spin_states1 = np.delete(self.state, [0])
        spin_states2 = np.delete(self.state, [self.length-1])
        self.energy = np.sum(self.bonds*spin_states1*spin_states2)
        return self
    
    #This finds the strongest bond
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
        self.local_energy = self.mega_bond*self.state[self.mega_index]*self.right_state # finds the energy that we will remove from the total energy during the transformation
        return self
    
    
    #This is the RG process
    def elimination_transformation(self):
        
        self.energy_prime = -(0.75)*self.mega_bond - ((0.1875)/self.mega_bond)*((self.left_mini_bond**2) + (self.right_mini_bond**2)) # find the energy contribution
        try:
            self.bond_prime = (self.left_mini_bond*self.right_mini_bond)/(2*self.mega_bond) # find the strength of the new bond that will exist after we remove the spins
        except IndexError:
            self.bond_prime = 0
        self.new_local_energy = self.energy_prime + (self.bond_prime*self.state[self.mega_index-1]*self.state[self.mega_index]) # finds the energy that we will add to the total energy after we remove the spins/bonds that existed
        
        self.state = np.delete(self.state, [self.mega_index, self.mega_index+1]) # the new chain after removing the spins sharing the strongest bond
        self.bonds[self.mega_index]=self.bond_prime # replacing the strongest bond with the bond prime
        self.bonds = np.delete(self.bonds, [self.mega_index-1, self.mega_index+1]) # removing the bonds next to the bond prime to get the proper new chain
        
        self.energy = self.energy - self.local_energy + self.new_local_energy # calculates the new energy of the chain by removing the previous contribution of the strongest bond and adding the new contribution of the newly weak bond in the same spot
        self.strong_bond() # refreshes the strongest bond
        self.mean_bonds()
        return self
        
