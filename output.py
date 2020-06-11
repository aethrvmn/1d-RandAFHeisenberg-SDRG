import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from spin_chains.zero_temp import ZT_Chain
from spin_chains.non_zero_temp import NZT_Chain

system2 = NZT_Chain(10000, 0.25, 0.1)

print(system2.logmega)

for i in tqdm(range(10)):
        system2.nzt_elimination()
print(system2.logmega)

        
