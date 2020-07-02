import matplotlib.pyplot as plt
import numpy as np
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(99,1,0.5)
matrix2 = NZT_Random_Spin(99,1,0.5, 1)

print(matrix1.left_index, matrix1.right_index )
matrix1.renormalization()
# # plt.matshow(matrix2.bond_matrix)
#
duration = time.time() - start_time
#print("---%s seconds---" %duration, matrix1.system_energy)
plt.matshow(matrix1.bond_matrix)
print(matrix1.left_index, matrix1.right_index )
plt.show()

