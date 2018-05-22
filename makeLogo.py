import gprpy as gp
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

filename = 'exampledata/SnS/ComOffs/XLINE00.DT1'

proj = gp.gprpy2d(filename)

maxpoint=100;

x=proj.twtt[0:maxpoint]

y=proj.data[0:maxpoint,10]

# Snake Body
lw=5
plt.plot(x,y,'k',linewidth=lw,solid_capstyle='round')

# Snake head
wid=2500
len1=-30
len2=10
xshift=0
Path = mpath.Path
path_data = [
    (Path.MOVETO, [xshift,wid]),
    (Path.CURVE3, [len1+xshift,0]),
    (Path.LINETO, [xshift,-wid]),
    (Path.CURVE3, [len2+xshift,0]),
    (Path.LINETO, [xshift,wid]),
    (Path.CLOSEPOLY, [xshift,wid])]
codes, verts = zip(*path_data)
path = mpath.Path(verts, codes)
patch = mpatches.PathPatch(path)
patch.set_facecolor('black')
plt.gca().add_patch(patch)

# Eyes
eye1 = mpatches.Ellipse([-2,1000], 3, 1000)
eye2 = mpatches.Ellipse([-2,-1000], 3, 1000)
eye1.set_facecolor('white')
eye2.set_facecolor('white')
plt.gca().add_patch(eye1)
plt.gca().add_patch(eye2)

# Tongue
x, y = np.array([[-10, -18, -20], [0.0, 0.0, 600]])
line1 = mlines.Line2D(x, y, lw=2, color='black')
x, y = np.array([[-10, -18, -20], [0.0, 0.0, -600]])
line2 = mlines.Line2D(x, y, lw=2, color='black')
plt.gca().add_line(line1)
plt.gca().add_line(line2)

plt.gca().set_xlim([-40,100])
#plt.gca().set_xlim([100,-40])
plt.gca().set_ylim([-20000,20000])


font = {'family': 'Verdana',
        'color':  'black',
        'weight': 'bold',
        'style': 'italic',
        'size': 35,
        }
plt.gca().text(35,-10000,'GPRPy',fontdict=font)

#plt.gca().axis('off')
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)

plt.savefig('GPRPy_Logo.pdf',format='pdf')
