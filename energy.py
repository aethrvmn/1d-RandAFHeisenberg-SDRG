import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin
from tqdm import tqdm

iterations=10

for j in tqdm(range(iterations)):
    matrix = ZT_Random_Spin(100,1)
    en = []
    i=0
    while matrix.end_rg == 0:
        matrix.renormalization()
        en = np.append(en, np.abs(matrix.system_energy))
        #print(matrix.system_energy)
    plt.plot(en)
plt.savefig('Figures/energy.png')
plt.show()