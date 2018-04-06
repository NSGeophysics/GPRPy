import gprpy as gp
import matplotlib.pyplot as plt
import numpy as np

filename = 'exampledata/SnS/ComOffs/XLINE00.DT1'
proj = gp.gprpy2d(filename)

#proj.setRange([0,32])

#proj.showHistory()

#print(len(proj.twtt))
#print(len(proj.profilePos))
#print(np.shape(proj.data))

#proj.writeHistory("histtest.py")

#proj.showTWTT(timelim=[0,200],profilelim=[50,200])
#ax = plt.ylim([0,200])
#ax.invert_yaxis()

#plt.imshow(proj.data)
#plt.show()



plt.show()
