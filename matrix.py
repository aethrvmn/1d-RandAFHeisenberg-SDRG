import matplotlib.pyplot as plt
from matplotlib import colors
import time
from random_singlets.NZT_RSS import NZT_Random_Spin
from random_singlets.ZT_RSS import ZT_Random_Spin

start_time = time.time()

matrix1 = ZT_Random_Spin(90,1,0.0001)
matrix2 = NZT_Random_Spin(99,1,0.5, 1)

cmap = colors.ListedColormap(['black','white','green'])
bounds=[-matrix1.ceiling, -0.0000000000001,0.0000000000001, matrix1.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)

while matrix1.max_bond > matrix1.floor and matrix1.end_rg == 0:
    matrix1.renormalization()
    img = plt.imshow(matrix1.bond_matrix,interpolation='nearest', cmap=cmap, norm=norm)
    plt.pause(0.01)
duration = time.time() - start_time #if anyone wants to flex their machine they should uncomment the print at the bottom.
plt.show()
#print("---%s seconds---" %duration) 

