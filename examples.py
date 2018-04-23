import gprpy as gp
import matplotlib.pyplot as plt
import numpy as np

filename = 'exampledata/SnS/ComOffs/XLINE00.DT1'
#filename = 'dewowed.gpr'
proj = gp.gprpy2d(filename)

proj.tpowGain(1.5)

#proj.showProfile()

#plt.show()

proj.printProfile('test1.pdf')

#proj.setVelocity(0.1)

#proj.printProfile('test2')


#proj.remMeanTrace(1000000)

#proj.dewow(100000000)

#proj.printProfile('testfig.pdf',timelim=[0,700])

#proj.writeHistory("test1.py")

#proj.undo()
#proj.undo()

#proj.writeHistory("test2.py")

#proj.timeZeroAdjust()
#proj.save('testsave')
#proj2 = gp.gprpy2d('./testsave.gpr')

#proj.toDepth(0.1)

#proj.showTWTT(timelim=[0,700])

#proj.printTWTT('testfig',timelim=[0,700])
#proj.prepTWTTfig(timelim=[0,700])
#plt.show()

#proj.dewow(100000000000000)
#proj.writeHistory("histtest.py")

#proj.timeZeroAdjust()
#proj.writeHistory("histtest1.py")
#proj.undo()
#proj2.writeHistory("histtest2.py")


#proj2.showTWTT()
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



#plt.show()
