import numpy as np

class spin_chain:

    '''
    1.0:
     This class creates instances of a spin chain, with a specific energy. For now generates the ground state of said system.
    When an instance is  created, the energy of the system is calculated automatically using the calculate_energy method, which is
    the typical way of calculating the energy with the regular hamiltonian, without spin elimination sheananigans (need to fix that?).
    The max_K method finds the maximum K and its index location and returns them as a tuple. To initialize the class you need to specify
    length (of lattice) and the probability cuttoff J for initialization of the K's. Also a state could be initialized if you want with
    different starting state in a newer version. Last and certainly not least the spin_elimination method... It is not perfect, and im
    obviously not sure if it does exactly what you want but i think it does... Also the case for the max K been found at the edges must
    be included. Now it does nothing (identity transform). For the regular case, the state from length goes to length - 2,
    (2 spins removed) as well as the K's change 3 die but one is added so you effectively lose 2 bonds, 1 from elimination, the
    other from the merging... (Hopefully correctly). The energy is also updated with hopefully the correct way ?! So when you want
    to do the spin elimination transform you, having initialized the class in an instance, namely idk bob = spin_chain(10, 1)
    you just do, bob.spin_elimination() and bob is updated !!! So you keep on doing that shit. Also every time you can access them
    K values, as: K_now = bob.K . Hopefully you get what i mean. Cheers!
    2.0:
     In this version i added the cases when K_max is located in the boundary. Not sure if correct but it works, and works in a manner
    that was anticipated, similar with the regular one.
    '''
    version = "2.0"

    def calculate_energy(self):
        #calculates energy of current state
        state = self.state

        s1 = np.delete(state, [0])
        s2 = np.delete(state, [self.length-1])

        return np.sum(s1*s2*self.K)

    def __init__(self, length, J):
        #initializes the groundd state using the specified length, as well as the random K's between [0, J)
        self.length = int(length)
        self.state = np.resize(np.array([0.5, -0.5]), int(length))
        self.K = J*np.random.rand(int(length) - 1)
        self.energy = self.calculate_energy()

    def max_K(self):

        K_max = np.amax(self.K)
        max_index = np.argmax(self.K)

        return K_max, max_index

    def spin_elimination(self):
        '''
        This function does the elimination transformation on the state and the energy as well as the K's
        '''
        K_old = self.K
        #max K and its index is found
        K_max, max_index = self.max_K()
        #i didn't know what happens if the transformation is on the boundaries :(
        if max_index == 0:
            K1 = 0
            K2 = self.K[max_index + 1]

            K_prime = K1*K2/(2*K_max)
            #updates all parameters that change!
            self.energy -= (K_max*(self.state[max_index]*self.state[max_index+1] + (3/4))) + (3*((K1**2)+(K2**2))/(16*K_max))

            K_old[max_index] = K_prime
            K_new = np.delete(K_old, [max_index, max_index+1])
            self.length = K_new.shape[0] + 1
            self.state = np.resize(np.array([0.5, -0.5]), self.length)
            self.K = K_new
        elif max_index == K_old.shape[0] - 1:
            K1 = self.K[max_index - 1]
            K2 = 0

            K_prime = K1*K2/(2*K_max)
            #updates all parameters that change!
            self.energy -= (K_max*(self.state[max_index]*self.state[max_index-1] + (3/4))) + (3*((K1**2)+(K2**2))/(16*K_max))

            K_old[max_index] = K_prime
            K_new = np.delete(K_old, [max_index-1, max_index])
            self.length = K_new.shape[0] + 1
            self.state = np.resize(np.array([0.5, -0.5]), self.length)
            self.K = K_new
        else:
            #the regular transformation
            K1 = self.K[max_index - 1]
            K2 = self.K[max_index + 1]

            K_prime = K1*K2/(2*K_max)
            #updates all parameters that change!
            self.energy -= (K_max*(self.state[max_index]*self.state[max_index+1] + (3/4))) + (3*((K1**2)+(K2**2))/(16*K_max))

            K_old[max_index] = K_prime
            K_new = np.delete(K_old, [max_index-1, max_index+1])
            self.length = K_new.shape[0] + 1
            self.state = np.resize(np.array([0.5, -0.5]), self.length)
            self.K = K_new
