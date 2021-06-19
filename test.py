import numpy as np
import matplotlib.pyplot as plt
from chain import Random_chain
from tqdm import tqdm

length = 1000
matrix = Random_chain(length,1)
eta = []

for i in tqdm(range(int(length/2))):
    eta = np.append(eta, matrix.eta)
    matrix.renormalization()

eta = np.sort(eta)
range=np.arange(len(eta))
sum = []

for i in tqdm(range):
    sum = np.append(sum, np.sum(eta[i:]))

plt.figure(1)
plt.plot(np.sort(sum/sum.max()), linewidth=1)
plt.plot(1-np.exp(-range), '--r', linewidth=1)
plt.ylim(0.9, 1.01)
plt.xlim(420, 500)
plt.yscale('log')
plt.yticks([1,1])
plt.ylabel('Value of $\eta$ near the random singlet phase (rescaled)')
plt.xlabel('Iteration step')
plt.savefig('Figures/Fixed/near_singlet1000.png')

plt.figure(2)
plt.plot(np.sort(sum/sum.max()), linewidth=1)
plt.plot(1-np.exp(-range), '--r', linewidth=1)
plt.ylabel('Value of $\eta$ during the SDRG (rescaled)')
plt.xlabel('Iteration step')
plt.yscale('log')
plt.savefig('Figures/Fixed/eta1000.png')