import numpy as np
import matplotlib.pyplot as plt
from chain import Random_chain
import seaborn as sns
import statsmodels.api as sm
from tqdm import tqdm

length=1000
val=np.zeros(length)
matrix = Random_chain(length,1,0)

leg=[]

plt.figure(1)
for i in tqdm(range(int(length/2))):
    if i == 0 or i==10 or i%100==0:
        plt.figure(1)
        sns.kdeplot(matrix.bond_matrix[matrix.bond_matrix != 0], bw_method=.5)
        leg = np.append(leg, 'Distribution for the '+str(i)+'th Iteration')
    matrix.renormalization()
    
plt.figure(2)
dens = sm.nonparametric.KDEUnivariate(matrix.bond_matrix[matrix.bond_matrix != 0])
dens.fit(bw=.5)
x =np.linspace(0,1,100)
y = dens.evaluate(x)

pp = np.sort(-np.log(np.abs(matrix.bond_matrix[matrix.bond_matrix != 0])))

step = 10
parts = int(len(pp)/step)
p=np.linspace(0,1,parts)
ppp = []
for i in tqdm(np.arange(0, len(pp), step)):
    ppp = np.append(ppp, np.sum(pp[i:]))
    
plt.figure(1)
plt.xlim(0, 1)
plt.ylabel('Distribution of bonds')
plt.xlabel('J')
plt.yticks([])
plt.legend(leg)
plt.savefig('Figures/distribution.png')

plt.figure(2)
plt.yticks([])
plt.ylabel('Distribution of bonds')
plt.plot(x,y)
plt.plot(x, np.exp(-x),'--')
plt.xlabel('J')
plt.savefig('Figures/singlet_dist.png')

plt.figure(3)
plt.plot(p,ppp)
plt.plot(p,np.exp(-p), '--')
plt.yscale('log')
plt.savefig('Figures/maybe(?).png')

plt.show()
