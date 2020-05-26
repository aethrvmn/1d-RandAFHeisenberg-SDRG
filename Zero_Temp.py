import numpy as np
import matplotlib.pyplot as plt
from classes.spin_chain_vasil import Chain
from classes.spin_chain_tolis import spin_chain

system1 = Chain(100,1)

for i in np.arange(100):
    system1.elimination_transformation()
    print(len(system1.bonds))
