import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from random_singlets.ZT_RSS import ZT_Random_Spin
import time

start_time = time.time()

N_systems = 1 #the amount of different systems we want to check
lattice_size = 30 # the number of spin particles
ceiling = 1 # the ceiling for the possible strength of the bonds
floor = 0.01 # the floor to stop the RG process at low energies
transformation_iterations = 300 # how many transformations should be made in each system
range_iterations = np.array(range(transformation_iterations))

plt.style.use('fivethirtyeight')


values = np.zeros(transformation_iterations)
log_values = np.zeros(transformation_iterations)
cov_values = np.zeros(transformation_iterations)

for i in tqdm(range(N_systems)):
    system = ZT_Random_Spin(lattice_size, ceiling, floor) #this creates the spin chain based on the chain class
    j=0
    while system.max_bond > system.floor:
        system.renormalization()
        values[j] = system.logmax
        j+=j
    log_values += values
    cov_values += values**2

log_values /= N_systems
cov_values /= N_systems
cov_values -= log_values**2
cov_values = np.sqrt(cov_values)

plt.errorbar(range_iterations, log_values, yerr = cov_values)
plt.plot(log_values)
plt.legend(loc = "upper right")
plt.ylabel("$\Gamma = -\ln(\Omega)$")
plt.xlabel("RG iterations")
plt.title(f"Change in $\Gamma$ for the Zero Temperature Case, N_systems = {N_systems}")
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
