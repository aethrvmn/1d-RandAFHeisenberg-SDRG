import matplotlib.pyplot as plt
import numpy as np
from random_chain import Random_chain
from scipy.stats import expon
from tqdm import tqdm
import seaborn as sns
import multiprocessing as mp

num_simulations = 100
num_cores = mp.cpu_count() // 2  # Use half of the available cores

def run_simulation(sim):
    length = 1000
    n_steps = int(length / 2)
    chain = Random_chain(length, resc = False)

    for _ in tqdm(range(n_steps)):
        chain.renormalization()

    final_bond_distribution = np.extract(chain.bond_matrix != 0, chain.bond_matrix)
    negative_bonds = [x for x in final_bond_distribution if x < 0]
    abs_negative_bonds = negative_bonds- min(negative_bonds)

    return abs_negative_bonds

if __name__ == '__main__':
    with mp.Pool(num_cores) as pool:
        bond_distributions = list(tqdm(pool.imap(run_simulation, range(num_simulations)), total=num_simulations))

    average_distribution = np.concatenate(bond_distributions)
    scale = np.mean(average_distribution)
    x = np.linspace(0, max(average_distribution), num=200)

    plt.figure(1)
    plt.hist(average_distribution, label="Average absolute values of negative bonds", density = True, bins = 'auto')
    plt.plot(x, expon.pdf(x, scale=scale), label="Exponential distribution", linestyle="--")
    plt.xlabel("Bond strength")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Average bond distribution")
    plt.savefig("Figures/average_bond_distribution_histogram.png")

    plt.figure(2)
    sns.kdeplot(average_distribution, label="Average absolute values of negative bonds", fill = True)
    plt.plot(x, expon.pdf(x, scale=scale), label="Exponential distribution", linestyle="--")
    plt.xlabel("Bond strength")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Average bond distribution")
    plt.savefig("Figures/average_bond_distribution_kdeplot.png")
    plt.show()
