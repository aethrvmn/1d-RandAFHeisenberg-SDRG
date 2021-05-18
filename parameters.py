import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin
from tqdm import tqdm
import seaborn as sns

length=100
iterations=20
genen=np.zeros(int(length/2))
genmean=np.zeros(int(length/2))
gendist=np.zeros(int(length/2))
for j in tqdm(range(iterations)):
    matrix = ZT_Random_Spin(length,1,0)
    en = []
    mean= []
    dist = []
    while matrix.end_rg == 0:
        matrix.renormalization()
        en = np.append(en, np.abs(matrix.system_energy))
        mean = np.append(mean, matrix.mean)
        dist = np.append(dist, matrix.distance)
    plt.figure(1)
    plt.plot(en, 'b')
    plt.figure(2)
    plt.plot(mean, 'b')
    plt.figure(3)
    #plt.hist(dist, bins=int(length/2))
    genen+=en
    genmean+=mean
    gendist+=dist
genen=genen/iterations
genmean = genmean/iterations
gendist = gendist/iterations

print(genen[int(length/2)-1])
plt.figure(1)
plt.plot(genen, '--r', linewidth=4)
plt.ylabel('System energy')
plt.xlabel('Iteration step')
plt.savefig('Figures/energy.png')

plt.figure(2)
plt.plot(genmean, '--r', linewidth=4)
plt.ylabel('Average bond strength')
plt.xlabel('Iteration step')
plt.savefig('Figures/meanstr.png')

plt.figure(3)
sns.kdeplot(gendist, bw_method=.3)
plt.xlim(0)
plt.ylabel('# of bonds')
plt.yticks([])
plt.xlabel('Distance between the two spins')
plt.savefig('Figures/distance.png')

plt.show()