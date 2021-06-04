import numpy as np
import matplotlib.pyplot as plt
from chain import Random_chain
from tqdm import tqdm
import seaborn as sns

length=1000
iterations=100

genen=np.zeros(int(length/2))
gendist=[]
genscale=np.zeros(int(length/2))

for j in tqdm(range(iterations)):
    matrix = Random_chain(length,1,0)
    en = []
    dist = []
    scale= []
    while matrix.end_rg == 0:
        matrix.parameters()
        en = np.append(en, np.abs(matrix.system_energy))
        dist = np.append(dist, matrix.distance)
        scale = np.append(scale, -np.log(matrix.max_bond))
        matrix.renormalization()
    plt.figure(1)
    plt.plot(en, 'b')
    plt.figure(2)
    plt.figure(3)
    plt.plot(scale, 'b')
    genen+=en
    gendist = np.append(gendist, dist)
    genscale+=scale

genen=genen/iterations
gendist = gendist/iterations
genscale = genscale/iterations

plt.figure(1)
plt.plot(genen, '--r', linewidth=2)
plt.ylabel('System energy')
plt.xlabel('Iteration step')
plt.yscale('log')
plt.savefig('Figures/energy.png')

plt.figure(2)
plt.hist(gendist, bins=int(length/2), histtype='step', log=True)
plt.ylabel('Number of bonds')
plt.xlabel('Distance between the two spins')
plt.savefig('Figures/distance.png')

plt.figure(3)
plt.plot(genscale, '--r', linewidth=2)
plt.ylabel('$-\ln(\Omega)$')
plt.xlabel('Iteration step')
plt.savefig('Figures/energyscale.png')
plt.yscale('log')
plt.xscale('log')

plt.show()
