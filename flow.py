import numpy as np
import matplotlib.pyplot as plt
from chain import Random_chain
from tqdm import tqdm

plt.ion

length = 1000
matrix = Random_chain(length)
eta = []

for i in tqdm(range(int(length/2))):
    eta = np.append(eta, matrix.eta)
    matrix.renormalization()

eta = np.sort(eta)
range=np.arange(len(eta))
sum = []

for i in range:
    sum = np.append(sum, np.sum(eta[i:]))

plt.figure(1)
plt.plot(np.sort(sum/sum.max()), linewidth=1)
plt.plot(1-np.exp(-range), '--r', linewidth=1)
plt.ylabel('Value of $\eta$ during the SDRG (rescaled)')
plt.xlabel('Iteration step')
plt.yscale('log')
plt.savefig('Figures/Fixed/eta_distribution.png')

plt.show()