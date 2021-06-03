import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
from chain import Random_chain

length = 100
matrix = Random_chain(length,1,0)

fig = plt.figure()

ax = plt.gca()
ax.set_xticks(np.arange(-0.5, length+1, 1))
ax.set_yticks(np.arange(-0.5, length+1, 1))
ax.set_xticklabels('')
ax.set_yticklabels('')
ax.grid(color='k', linestyle='-', linewidth=0.05)

cmap = colors.ListedColormap(['black','white','yellow'])
bounds=[-matrix.ceiling, -0.0000000000001,0.0000000000001, matrix.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)
ims = []
img = plt.imshow(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
img.figure.savefig('Figures/Matrix/startmatrix.png', pad_inches=0, transparent=True)

i=0
while matrix.end_rg == 0:
    matrix.renormalization()
    if i == 1 or i == 5 or i == 10 or i ==13:
        img.figure.savefig('Figures/Matrix/midrg'+str(i)+'.png', pad_inches=0, transparent=True)
    img = plt.imshow(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
    ims.append([img])
    i+=1

ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=10000)
ani.save('Figures/rg.gif')

img.figure.savefig('Figures/Matrix/singlet-phase.png', pad_inches=0, transparent=True)
plt.show()
