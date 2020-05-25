import numpy as np
import matplotlib.pyplot as plt
from classes.spin_chain_vasil import Chain
from classes.spin_chain_tolis import spin_chain

system = Chain(100,10)

system.strong_bond().sys_energy()

print(system.mega_bond, system.mega_index, system.energy)
