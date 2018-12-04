import sys
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import simpledialog as sd
from tkinter import messagebox as mesbox
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import gprpy as gp
import numpy as np
#import gprpy.toolbox.splash as splash
import os
import Pmw
import scipy.interpolate as interp

from gprpy.gprpyGUI import GPRPyApp
from gprpy.gprpyCWGUI import GPRPyCWApp



if len(sys.argv) < 2:
    mode = input("Profile [p] or Common Midpoint / WARR [c]?  ")
else:
    mode = sys.argv[1][0]



if mode == 'p':
    #colsp=2
    rightcol=9
    #halfwid=6
    figrowsp=19+1
    #figcolsp=9
    
    root = tk.Tk()
    
    for col in range(rightcol):
        root.columnconfigure(col, weight=1)
    for row in range(figrowsp):    
        root.rowconfigure(row, weight=1)
            
    app = GPRPyApp(root)

    root.mainloop()

elif mode == 'c' or mode == 'w':
    #colsp=2
    rightcol=10
    #halfwid=4
    figrowsp=15+1
    #tagx=10 # 2
    #tagy=5 # -3
    
    root = tk.Tk()

    for col in range(rightcol):
        root.columnconfigure(col, weight=1)
    for row in range(figrowsp):    
        root.rowconfigure(row, weight=1)

    app = GPRPyCWApp(root)

    root.mainloop()

else:
    print("You need to chose either profile [p] or CMP/WARR [c] mode.")
    
    

