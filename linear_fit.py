import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from classes.zero_temp_spin_chain import Chain
from scipy.optimize import curve_fit

N_systems = 1000 #the amount of different systems we want to check
lattice_size = 1001 # the number of spin particles
ceiling = 100 # the ceiling for the possible strength of the bonds
transformation_iterations = 100 # how many transformations should be made in each system

# The linear fit function
def linear(x, a, b):
    y = a*x + b
    return y



for i in range(N_systems)
    system = Chain(lattice_size, ceiling) #this creates the spin chain based on the chain class
    for j in range(transformation_iterations):
        system.elimination_transformation()
