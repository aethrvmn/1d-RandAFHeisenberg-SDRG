import numpy as np


'''
This is the code I have writen to define the Random Antiferromagnetic Heisenberg Spin Chain and the RG transformation (elimination transformation) for a zero Temprature system.
'''
class Chain:
    
    version='Alpha v0.1'
    
    #Initial call, starts the project by taking in the length and the maximum energy ceiling, and creating a random spin chain based on those.
    def __init__(self, length, ceiling):
        self.length = int(length)
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.bonds = ceiling*np.random.rand(int(length))
        
        
    #This calculates the system energy    
    def sys_energy(self):
        #as of right now does not take into account the last spin as it IndexErrors
        for i in np.arange(self.length):
            try:
                energy =+ self.bonds[i]*self.state[i]+self.state[i+1]
                continue
            except IndexError:
                break
        return energy
    
    
    #This finds the strongest bond
    def biggest_bond(self):
        mega_bond = np.amax(self.bonds)
        mega_index = np.argmax(self.bonds)
        return mega_bond, mega_index
    
    
    #This is the RG process, does not account for the edge spins
    def elimination_transformation(self):
        
