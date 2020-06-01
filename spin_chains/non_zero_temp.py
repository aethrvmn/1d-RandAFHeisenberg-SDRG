import numpy as np


'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain and the RG transformation (elimination transformation) for a non-zero temprature system. Currently there is an issue with what happens when the size of the chain falls to zero, and shows an error. I also need to change the distribution method for the bonds so it's an actual random distribution to properly study the RG flow.
'''
class Non_Zero_Chain:
    
    version = 'v0.1'
    
    
    # Initial call, starts the project by taking in the length and the maximum energy ceiling, and creating a random spin chain based on those.
    def __init__(self, length, ceiling, temprature):#, floor):
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.bonds = ceiling*np.random.rand(len(self.state) - 1)
        self.length = len(self.state)
        #self.floor = floor
        self.sys_energy() # compute the energy
        self.strong_bond() # find the strongest bond
        self.mean_bonds() # compute the mean 
        
    # Calculates the mean value of the bonds    
    def mean_bonds(self):
        self.mean = np.mean(self.bonds)
        
    #Instead of the ground energy, since we have a finite temperature, we calculate the free energy
