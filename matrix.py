import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(9,1,0.5)
matrix2 = NZT_Random_Spin(99,1,0.5, 1)

cmap = colors.ListedColormap(['black','white','blue'])
bounds=[-matrix1.ceiling, -0.0000000001,0.0000000001, matrix1.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)


#matrix1.renormalization()
# # plt.matshow(matrix2.bond_matrix)

duration = time.time() - start_time
#print("---%s seconds---" %duration, matrix1.system_energy)
img = plt.imshow(matrix1.bond_matrix,interpolation='nearest', cmap=cmap, norm=norm)
#print(matrix1.left_index, '\n',matrix1.right_index )
#plt.show()

