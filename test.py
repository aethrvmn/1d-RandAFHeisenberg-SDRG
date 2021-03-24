import numpy as np
import matplotlib.pyplot as plt
from random_singlets.ZT_RSS import ZT_Random_Spin

matrix = ZT_Random_Spin(30,1)

values = np.zeros(15)
itrange=np.arange(15)
i=0

while matrix.end_rg == 0:
    matrix.renormalization()
    matrix.rgflow()
    #print(-np.log(np.abs(np.prod(matrix.flow))))
    values[i]=-np.log(np.abs(np.prod(matrix.flow)))
    i+=1
plt.plot(itrange, values, 'b')
#plt.fill_between(itrange, values, color='blue', alpha=0.3)


#print(len(matrix.flow)-len(values))
plt.show()
#while matrix.end_rg == 0:
#    matrix.renormalization()
    
