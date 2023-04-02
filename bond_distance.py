import numpy as np
import matplotlib.pyplot as plt
from random_chain import Random_chain
from tqdm import tqdm
import seaborn as sns
import multiprocessing as mp

length = 1000
num_simulations = 100
num_cores = mp.cpu_count() // 2  # Use half of the available cores

# Initialize an empty list to store the bond distances from each simulation
all_distances = []

def run_simulation(sim):
    matrix = Random_chain(length)
    distances = []

    for _ in tqdm(range(int(length / 2))):
        matrix.renormalization()

        # Calculate the distance between the two spins
        distance = matrix.right_index[1] - matrix.left_index[1] - 1
        distances.append(distance)

    return distances

if __name__ == '__main__':
    with mp.Pool(num_cores) as pool:
        distance_lists = list(tqdm(pool.imap(run_simulation, range(num_simulations)), total=num_simulations))

    # Flatten the list of lists into a single list
    all_distances = [item for sublist in distance_lists for item in sublist]

    # Plot the kdeplot for the bond distances
    plt.hist(all_distances, label="Distances between bonds", density = True, bins = 'auto')
    plt.yscale('log')
    plt.ylabel('Density')
    plt.xlabel('Number of In-Between Spins')
    plt.legend()
    plt.savefig('Figures/distance.png')

    plt.show()
