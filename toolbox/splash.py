import gprpy as gp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import os

def showSplash(a,dir_path):
    filename=os.path.join(dir_path,'exampledata','SnS','ComOffs','XLINE00.DT1')
    snakeGPR = gp.gprpy2d(filename)
    maxpoint=100;
    x=snakeGPR.twtt[0:maxpoint]
    y=snakeGPR.data[0:maxpoint,10]
    # Snake body
    lw=5
    a.plot(x,y,'k',linewidth=lw,solid_capstyle='round')
    # Snake head
    Path = mpath.Path
    xshift=0
    path_data = [
        (Path.MOVETO, [xshift,2500]),
        (Path.CURVE3, [-30+xshift,0]),
        (Path.LINETO, [xshift,-2500]),
        (Path.CURVE3, [10+xshift,0]),
        (Path.LINETO, [xshift,2500]),
        (Path.CLOSEPOLY, [xshift,2500])]
    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path)
    patch.set_facecolor('black')
    a.add_patch(patch)
    # Eyes
    eye1 = mpatches.Ellipse([-2,1000], 3, 1000)
    eye2 = mpatches.Ellipse([-2,-1000], 3, 1000)
    eye1.set_facecolor('white')
    eye2.set_facecolor('white')
    a.add_patch(eye1)
    a.add_patch(eye2)
    # Tongue
    x, y = np.array([[-10, -18, -20], [0.0, 0.0, 600]])
    line1 = mlines.Line2D(x, y, lw=2, color='black')
    x, y = np.array([[-10, -18, -20], [0.0, 0.0, -600]])
    line2 = mlines.Line2D(x, y, lw=2, color='black')
    a.add_line(line1)
    a.add_line(line2)
    # Axis setup
    a.set_xlim([-40,100])
    a.set_ylim([-20000,20000])
    #a.axis('off')
    # Text
    font = {'family': 'Verdana',
        'color':  'black',
        'weight': 'bold',
        'style': 'italic',
        'size': 35,
        }
    a.text(35,-10000,'GPRPy',fontdict=font)
