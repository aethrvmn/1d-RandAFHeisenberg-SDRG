import numpy as np
import matplotlib.pyplot as plt
from classes.spin_chain_vasil import Chain

system1 = Chain(10000,1)

for i in np.arange(1000):
    system1.elimination_transformation()
    plt.scatter(len(system1.bonds), system1.mean)
plt.show()
