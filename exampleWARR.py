import gprpy
mygpr = gprpy.gprpyCW()
mygpr.importdata('/Users/aplattner/Desktop/MySoftware/Git/GPRPy/exampledata/SnS/WARR/XLINE00.DT1',dtype='WARR')
mygpr.normalize()

vmin=0.05
vmax=0.35
vint=0.01

mygpr.linSemblance(vmin,vmax,vint)

import matplotlib.pyplot as plt
import numpy as np

v = np.arange(vmin,vmax+vint,vint)
plt.imshow(np.abs(np.asmatrix(mygpr.linSemb)).transpose(),extent=[vmin,vmax,np.max(mygpr.twtt),np.min(mygpr.twtt)],aspect='auto')
plt.colorbar()
plt.show()
