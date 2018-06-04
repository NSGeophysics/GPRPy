import gprpy
import numpy as np
import toolbox.gprpyTools as tools
#import pyximport
#pyximport.install(setup_args={"include_dirs":np.get_include()})
#import toolbox.gprpyTools2 as tools2

import matplotlib.pyplot as plt


mygpr = gprpy.gprpyCW()
mygpr.importdata('/Users/aplattner/Desktop/MySoftware/Git/GPRPy/exampledata/SnS/WARR/XLINE00.DT1',dtype='WARR')
#mygpr.normalize()

mygpr.setZeroTime(31)
#mygpr.showProfile()

#plt.imshow(mygpr.data)

fig = plt.figure()


vmin=0.03
vmax=0.35
vint=0.005
vVals = np.arange(vmin,vmax,vint)
tVals = np.arange(1.18,700.0,0.8)


# linSemb = tools.linSemblance_alt2(mygpr.data,mygpr.profilePos,mygpr.twtt,vVals,tVals,typefact=1.0)
# ax1 = fig.add_subplot(1,2,1)
# ax1.imshow(linSemb,extent=[vmin,vmax,max(tVals),min(tVals)],aspect='auto')

# linSembC = tools2.linSemblance(mygpr.data,mygpr.profilePos,mygpr.twtt,vVals,tVals,typefact=1.0)
# ax1 = fig.add_subplot(1,2,1)
# ax1.imshow(linSembC,extent=[vmin,vmax,max(tVals),min(tVals)],aspect='auto')

# linSemb2 = tools.linSemblance(mygpr.data,mygpr.profilePos,mygpr.twtt,vVals,tVals,typefact=1.0)
# ax2 = fig.add_subplot(1,2,2)
# ax2.imshow(linSemb2,extent=[vmin,vmax,max(tVals),min(tVals)],aspect='auto')



hypSemb = tools.hypSemblance(mygpr.data,mygpr.profilePos,mygpr.twtt,vVals,tVals,typefact=1.0)
#ax2 = fig.add_subplot(1,1,1)
plt.imshow(hypSemb,extent=[vmin,vmax,max(tVals),min(tVals)],aspect='auto')

plt.colorbar()
plt.show()
