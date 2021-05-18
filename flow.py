import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin
import seaborn as sns
import statsmodels.api as sm
from tqdm import tqdm

length=1000
val=np.zeros(length)
matrix = ZT_Random_Spin(length,1,0)

leg=[]

plt.figure(1)
for i in tqdm(range(int(length/2))):
    if i == 0 or i==10 or i==200 or i==int(length/2 -1):
        plt.figure(1)
        sns.kdeplot(matrix.bond_matrix[matrix.bond_matrix != 0], bw_method=.5)
        leg = np.append(leg, 'Distribution for the '+str(i)+'th Iteration')

    matrix.renormalization()
plt.figure(2)
dens = sm.nonparametric.KDEUnivariate(matrix.bond_matrix[matrix.bond_matrix != 0])
dens.fit(bw=.5)
x =np.linspace(0,1,100)
y = dens.evaluate(x)

plt.figure(1)
plt.xlim(0, 1)
plt.xlabel('')
plt.ylabel('')
plt.yticks([])
plt.legend(leg)
plt.savefig('Figures/distribution.png')

plt.figure(2)
plt.yticks([])
plt.plot(x,y)
plt.xlabel('')
plt.savefig('Figures/singlet_dist.png')

plt.show()