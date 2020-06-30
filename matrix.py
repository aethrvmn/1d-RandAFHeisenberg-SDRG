import matplotlib.pyplot as plt
import numpy as np
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix = ZT_Random_Spin(1000,100,0.5)

plt.matshow(matrix.bond_matrix)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
