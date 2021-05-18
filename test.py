import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin
from tqdm import tqdm

length=1000
matrix = ZT_Random_Spin(length,1,0)

for i in tqdm(range(int(length/2)-1)):
    matrix.renormalization()
dens = sm.nonparametric.KDEUnivariate(np.abs(matrix.bond_matrix[matrix.bond_matrix != 0]))
dens.fit()
x =np.linspace(0,1,100) #restrict range to (0,1)
y = dens.evaluate(x)
plt.plot(x,y, 'k')
plt.yticks([])
plt.xlabel('')
plt.savefig('Figures/test.png')