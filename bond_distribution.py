import matplotlib.pyplot as plt
import numpy as np
from random_chain import Random_chain
from scipy.stats import expon
from tqdm import tqdm
import seaborn as sns

num_simulations = 20

bond_distributions = []

length = 1000
n_steps = int(length / 2)

for sim in tqdm(range(num_simulations)):
    chain = Random_chain(length, resc = False)

    for _ in range(n_steps):
        chain.renormalization()

    # Calculate the final bond distribution
    final_bond_distribution = np.extract(chain.bond_matrix != 0, chain.bond_matrix)

    # Separate positive and negative bonds
    negative_bonds = [x for x in final_bond_distribution if x < 0]

    # Get the absolute values of the negative bonds
    abs_negative_bonds = negative_bonds - min(negative_bonds)

    # Calculate the scale parameter for the exponential distribution
    bond_distributions.append(abs_negative_bonds)

# Combine the bond distributions from each run and calculate the average distribution
average_distribution = np.concatenate(bond_distributions)

# Calculate the scale parameter for the exponential distribution
scale = np.mean(average_distribution)

# Generate x values for the exponential distribution
x = np.linspace(0, max(average_distribution), num=200)

# Plot the histogram for the average distribution
plt.hist(average_distribution, bins='auto', density=True, alpha=0.7, label="Average absolute values of negative bonds")

# Plot the exponential distribution
plt.plot(x, expon.pdf(x, scale=scale), label="Exponential distribution", linestyle="--")
plt.xlabel("Bond strength")
plt.ylabel("Density")
plt.legend()
plt.title("Average bond distribution")
plt.savefig("Figures/average_bond_distribution_histogram.png")
plt.show()

