import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from spin_chains.zero_temp import ZT_Chain
from spin_chains.non_zero_temp import NZT_Chain

#system2 = NZT_Chain(10000, 0.25, 0.1)
for j in range(10):
    system1 = ZT_Chain(100000, 0.25)
    system2 = NZT_Chain(100000, 0.25, 0.1)

    for i in tqdm(range(10000)):
        system1.zt_elimination()
        system2.nzt_elimination()
        
