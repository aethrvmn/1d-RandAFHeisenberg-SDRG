import numpy as np
import matplotlib.pyplot as plt
from random_chain import Random_chain
from tqdm import tqdm
import seaborn as sns

length = 1000
num_simulations = 20

# Initialize an empty list to store the bond distances from each simulation
all_distances = []

for sim in tqdm(range(num_simulations)):
    matrix = Random_chain(length)
    distances = []

    for _ in range(int(length / 2)):
        matrix.renormalization()

        # Calculate the distance between the two spins
        distance = matrix.right_index[1] - matrix.left_index[1] - 1
        distances.append(distance)

    # Combine the distances from all simulations
    all_distances.extend(distances)

# Plot the kdeplot for the bond distances
sns.kdeplot(all_distances, label="Distances between bonds", fill=True)
plt.yscale('log')
plt.ylabel('Density')
plt.xlabel('Number of In-Between Spins')
plt.legend()
plt.savefig('Figures/distance.png')

plt.show()

