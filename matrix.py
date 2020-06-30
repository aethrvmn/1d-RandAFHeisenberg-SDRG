import matplotlib.pyplot as plt
import numpy as np
from random_singlets.ZT_RSS import ZT_Random_Spin

matrix = ZT_Random_Spin(100,1,0.5)

plt.matshow(matrix.bond_matrix)
plt.show()
