import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
from random_singlets.ZT_RSS import ZT_Random_Spin

fig = plt.figure()

matrix = ZT_Random_Spin(30,1)

cmap = colors.ListedColormap(['black','white','yellow'])
bounds=[-matrix.ceiling, -0.0000000000001,0.0000000000001, matrix.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)
ims = []

while matrix.end_rg == 0:
    matrix.renormalization()
    img = plt.imshow(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
    ims.append([img])
    #plt.pause(1)

ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=10000)
ani.save('rg.mp4')
img.figure.savefig('singlet-phase.png')
plt.show()



