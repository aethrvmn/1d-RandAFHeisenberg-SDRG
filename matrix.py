import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
from random_singlets.ZT_RSS import ZT_Random_Spin

fig = plt.figure()

matrix = ZT_Random_Spin(50,1)

cmap = colors.ListedColormap(['black','white','yellow'])
bounds=[-matrix.ceiling, -0.0000000000001,0.0000000000001, matrix.ceiling]
norm = colors.BoundaryNorm(bounds, cmap.N)
ims = []
dist = []
i=0
while matrix.end_rg == 0:
    matrix.renormalization()
    #img = plt.hist(matrix.bond_matrix, interpolation='nearest', cmap=cmap, norm=norm, animated=True)
    #ims.append([img])
    dist = np.append(dist, matrix.distance)
    #plt.pause(1)

histo = plt.hist(dist)
#ani = animation.ArtistAnimation(fig, ims, interval=5000, blit=True, repeat_delay=10000)
#histo.savefig('distance.png', img)
#ani.save('woah.mp4')
plt.show()



