import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from random_singlets.ZT_RSS import ZT_Random_Spin
import time

#start_time = time.time()

#N_systems = 1000 #the amount of different systems we want to check
lattice_size = 90 # the number of spin particles
ceiling = 1 # the ceiling for the possible strength of the bonds
#floor = 0.000000001 # the floor to stop the RG process at low energies
#Num_iteration = 3000 # how many transformations should be made in each system


Num_models = 1000
Num_iteration = 40
range_iterations = np.array(range(Num_iteration))

values = np.zeros(Num_iteration)
log_values = np.zeros(Num_iteration)
cov_values = np.zeros(Num_iteration)

for i in tqdm(range(Num_models)):
    system = ZT_Random_Spin(lattice_size, ceiling)
    for j in range(Num_iteration):
        log_values[j] = system.logmax
        system.renormalization()
    values += log_values
    cov_values += log_values**2

values /= Num_models
cov_values /= Num_models
cov_values -= values**2
cov_values = np.sqrt(cov_values)

plt.errorbar(range_iterations, values, yerr = cov_values)
#plt.plot(values)
plt.legend(loc = "upper right")
plt.ylabel("$\Gamma = -\ln(\Omega)$")
plt.xlabel("RG iterations")
plt.title(f"Change in $\Gamma$ for the Zero Temperature Case, N_systems = {Num_models}")
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
