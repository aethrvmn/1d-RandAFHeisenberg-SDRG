import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin

matrix = ZT_Random_Spin(100,1)

en = []

while matrix.end_rg == 0:
    matrix.renormalization()
    en = np.append(en, -matrix.system_energy)

plt.plot(en)
plt.savefig('energy.png')
plt.show()
