import gprpy as gp
import matplotlib.pyplot as plt


filename = 'exampledata/SnS/WARR/XLINE00.DT1'
proj = gp.gprpy2d(filename)

plt.imshow(proj.data)
plt.show()

