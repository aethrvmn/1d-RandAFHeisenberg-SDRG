import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin

matrix = ZT_Random_Spin(100,1)

dist = []

while matrix.end_rg == 0:
    matrix.renormalization()
    dist = np.append(dist, matrix.distance)

plt.hist(dist)
plt.savefig('distance.png')
plt.show()
