import gprpy as gp
import matplotlib.pyplot as plt
import numpy as np

filename = 'exampledata/SnS/ComOffs/XLINE00.DT1'
proj = gp.gprpy2d(filename)

#proj.showTWTT(timelim=[0,200])
#plt.show()



proj.timeZeroAdjust()
proj.writeHistory("histtest1.py")
proj.undo()
proj.writeHistory("histtest2.py")


#proj.showTWTT(timelim=[0,200])
#plt.show()






#proj.setRange([0,32])

#proj.showHistory()

#print(len(proj.twtt))
#print(len(proj.profilePos))
#print(np.shape(proj.data))

#proj.writeHistory("histtest.py")


#ax = plt.ylim([0,200])
#ax.invert_yaxis()

#plt.imshow(proj.data)
#plt.show()



plt.show()
