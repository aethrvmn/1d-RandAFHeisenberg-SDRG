import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from spin_chains.zero_temp import ZT_Chain
from spin_chains.non_zero_temp import NZT_Chain

system2 = NZT_Chain(10000, 0.25, 0.001)
#for j in range(10):
    #system1 = ZT_Chain(100000 ,100)

    #for i in tqdm(range(10000)):
        #system1.elimination_transformation()

#print(len(system1.bonds))
system2.sys_free_energy()
print(system2.free_zero)
