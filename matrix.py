import matplotlib.pyplot as plt
import numpy as np
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(100,10,0.5)
matrix2 = NZT_Random_Spin(100,10,0.5, 1)


print(matrix1.logbonds)
# matrix1.renormalization()
# # plt.matshow(matrix1.bond_matrix)
# # plt.matshow(matrix2.bond_matrix)
# # plt.show()
#
duration = time.time() - start_time
print("---%s seconds---" %duration)
