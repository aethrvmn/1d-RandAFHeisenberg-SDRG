import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap, BoundaryNorm
from random_chain import Random_chain
from tqdm import tqdm

length = 30
matrix = Random_chain(length)

fig, ax = plt.subplots()

ax.set_xticks(np.arange(-0.5, length+1, 1))
ax.set_yticks(np.arange(-0.5, length+1, 1))
ax.set_xticklabels('')
ax.set_yticklabels('')
ax.grid(color='k', linestyle='-', linewidth=0.05)

cmap = ListedColormap(["black", "white", "yellow"])
norm = BoundaryNorm([-matrix.ceiling, -1e-10, 1e-10, matrix.ceiling], cmap.N)

ims = []
img = ax.imshow(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
fig.savefig('Figures/Matrix/startmatrix.png', pad_inches=0, transparent=True)

for i in tqdm(range(int(length/2))):
    matrix.renormalization()
    img = ax.imshow(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
    ims.append([img])

    if i == 1 or i == 5 or i == 10 or i == 13:
        fig.savefig('Figures/Matrix/midrg'+str(i)+'.png', pad_inches=0, transparent=True)

fig.savefig('Figures/Matrix/singlet-phase.png', pad_inches=0, transparent=True)

ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=10000)
ani.save('Figures/rg.gif')
plt.show()
