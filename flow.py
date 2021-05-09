import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin
from tqdm import tqdm

iterations=10

for j in tqdm(range(iterations)):
    matrix = ZT_Random_Spin(100,1)
    values = []
    i=0
    while matrix.end_rg == 0:
        matrix.renormalization()
        matrix.rgflow()
        values = np.append(values, np.sum(matrix.flow[matrix.flow > 0])/matrix.length)
        i+=1
    plt.plot(values)

plt.savefig('Figures/rgflow.png')
plt.show()