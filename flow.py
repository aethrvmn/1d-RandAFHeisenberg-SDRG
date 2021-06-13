import numpy as np
import matplotlib.pyplot as plt
from chain import Random_chain
import seaborn as sns
import statsmodels.api as sm
from tqdm import tqdm

def generator():
    while matrix.end_rg == 0:
        yield

length = 200
matrix = Random_chain(length,1)
eta = []
leg=[]
i=0

for _ in tqdm(generator()):
    if i == 0 or i==10 or i%40==0:
        plt.figure(1)
        sns.kdeplot(matrix.bond_matrix[matrix.bond_matrix != 0], bw_method=.5)
        leg = np.append(leg, 'Distribution for the '+str(i)+'th Iteration')
    eta = np.append(eta, matrix.eta)
    matrix.renormalization()
    i+=1

step = 1
parts = int(len(eta)/step)
sum = []

for i in tqdm(np.arange(0, len(eta), step)):
    val = np.sum(eta[i:])
    if val != 0:
        sum = np.append(sum, val)

rescaled_range=np.linspace(0,1,len(sum))

plt.figure(1)
plt.figure(1)
plt.xlim(0, 1)
plt.ylabel('Distribution of bonds')
plt.xlabel('Strength $J$ of the bonds')
plt.legend(leg)
plt.savefig('Figures/distribution.png')

plt.figure(2)

plt.plot(rescaled_range, sum/sum.max())
plt.plot(rescaled_range, np.exp(-5*rescaled_range), '--r')
plt.yscale('log')
plt.ylabel('Value of $\eta$ during the SDRG')
plt.xlabel('Iteration step (rescaled).')
plt.savefig('Figures/fixed_dist.png')

plt.show()
