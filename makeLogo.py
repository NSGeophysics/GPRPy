import gprpy as gp
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import matplotlib.image as im

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

plt.gca().set_xlim([-25,90])
#plt.gca().set_xlim([100,-40])
plt.gca().set_ylim([-28000,12000])


font = {'family': 'DejaVu Sans',
        'color':  'black',
        'weight': 'bold',
        'style': 'italic',
        'size': 35,
        }
plt.gca().text(35,-10000,'GPRPy',fontdict=font)

plt.gca().axis('off')
#plt.gca().get_xaxis().set_visible(False)
#plt.gca().get_yaxis().set_visible(False)

# add nsf logo
nsf =im.imread('toolbox/splashdat/NSF_4-Color_bitmap_Logo.png')
yanchor = -25000
yheight = 10000
xanchor = -20
#ratio = plt.gca().get_xlim/plt.gca().get_ylim
ratio = plt.gca().get_data_ratio()*1.36
xwidth = yheight/ratio
plt.gca().imshow(nsf, aspect='auto', extent=(xanchor, xanchor+xwidth, yanchor, yanchor+yheight))
font2 = {'family': 'DejaVu Sans',
        'color':  'black',
        'size': 9.5
        }
plt.gca().text(-20,-27000,'EAR-1550732',fontdict=font2)

plt.savefig('GPRPy_Logo.pdf',format='pdf',dpi=600)
#plt.savefig('NSF_Logo.pdf',format='pdf')
