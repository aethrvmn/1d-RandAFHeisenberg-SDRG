import numpy as np
import matplotlib.pyplot as plt

'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain and the RG transformation (elimination transformation) for a zero temprature system.
'''
class Chain:
    
    version='v1.1'
    

    #Initial call, starts the project by taking in the length and the maximum energy ceiling, and creating a random spin chain based on those.
    def __init__(self, length, ceiling):
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.bonds = ceiling*np.random.rand(len(self.state) - 1)
        self.length = len(self.state)
        self.sys_energy()
    
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
        try:
            self.left_mini_bond = self.bonds[self.mega_index-1] # the bond between the left spin and the one left to it
        except IndexError:
            self.left_mini_bond = 0
        try:
            self.right_mini_bond = self.bonds[self.mega_index+1] # the bond between the right spin and the one to the right of it
        except IndexError:
            self.right_mini_bond = 0
        self.local_energy = self.mega_bond*self.state[self.mega_index]*self.state[self.mega_index+1] # finds the energy that we will remove from the total energy during the transformation
        return self
    
    
    #This is the RG process
    def elimination_transformation(self):
        self.strong_bond() # compute the energy and find the strongest bond
        
        self.energy_prime = -(0.75)*self.mega_bond - ((0.1875)/self.mega_bond)*((self.left_mini_bond**2) + (self.right_mini_bond**2)) # find the energy contribution
        self.bond_prime = (self.left_mini_bond*self.right_mini_bond)/(2*self.mega_bond) # find the strength of the new bond that will exist after we remove the spins
        self.new_local_energy = self.energy_prime + (self.bond_prime*self.state[self.mega_index-1]*self.state[self.mega_index+2]) # finds the energy that we will add to the total energy after we remove the spins/bonds that existed
        
        self.energy = self.energy - self.local_energy + self.new_local_energy # calculates the new energy of the chain by removing the previous contribution of the strongest bond and adding the new contribution of the newly weak bond in the same spot

        
        self.state = np.delete(self.state, [self.mega_index, self.mega_index+1]) # the new chain after removing the spins sharing the strongest bond
        self.bonds[self.mega_index]=self.bond_prime # replacing the strongest bond with the bond prime
        self.bonds = np.delete(self.bonds, [self.mega_index-1, self.mega_index+1]) # removing the bonds next to the bond prime to get the proper new chain
        self.strong_bond() # refreshes the strongest bond
        return self
        
