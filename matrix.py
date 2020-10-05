import matplotlib.pyplot as plt
from matplotlib import colors
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(30,1,0.0000001)
matrix2 = NZT_Random_Spin(99,1,0.5, 1)

cmap = colors.ListedColormap(['black','white','green'])
bounds=[-matrix1.ceiling, -0.0000000000001,0.0000000000001, matrix1.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)

while matrix1.max_bond > matrix1.floor:
    matrix1.renormalization()
    #print(matrix1.left_index, matrix1.max_index, matrix1.right_index)
    print(matrix1.max_bond)
    #time.sleep(1)
    img = plt.imshow(matrix1.bond_matrix,interpolation='nearest', cmap=cmap, norm=norm)
    #plt.pause(0.0001)
plt.show()
duration = time.time() - start_time
print("---%s seconds---" %duration)

