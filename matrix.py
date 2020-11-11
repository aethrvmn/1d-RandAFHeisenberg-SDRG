#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from random_singlets.ZT_RSS import ZT_Random_Spin


matrix = ZT_Random_Spin(90,1)

cmap = colors.ListedColormap(['black','white','yellow'])
bounds=[-matrix.ceiling, -0.0000000000001,0.0000000000001, matrix.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)
#typ_length=np.array([])
#true_length=np.array([])
while matrix.end_rg == 0:
    matrix.renormalization()
    img = plt.imshow(matrix.bond_matrix,interpolation='nearest', cmap=cmap, norm=norm)
    plt.pause(0.01)
    #typ_length = np.append([matrix.logmax**2])
    #true_length = np.append([matrix.max_index[0] -matrix.max_index[1]])
plt.show()

