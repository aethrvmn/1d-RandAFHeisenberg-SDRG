import matplotlib.pyplot as plt
from matplotlib import colors
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(20,1,0.5)
matrix2 = NZT_Random_Spin(99,1,0.5, 1)

cmap = colors.ListedColormap(['black','white','yellow'])
bounds=[-matrix1.ceiling, -0.0000000000001,0.0000000000001, matrix1.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)

for i in range(10):
    matrix1.renormalization()
    print(matrix1.left_index, matrix1.max_index, matrix1.right_index)
    img = plt.imshow(matrix1.bond_matrix,interpolation='nearest', cmap=cmap, norm=norm)
    plt.show()
duration = time.time() - start_time
print("---%s seconds---" %duration)

