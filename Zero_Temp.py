import numpy as np
import matplotlib.pyplot as plt
from classes.spin_chain_vasil import Chain
from classes.spin_chain_tolis import spin_chain

system1 = Chain(100,1)

system1.strong_bond().sys_energy()
print(system1.energy)
system1.elimination_transformation()
print(system1.energy)


