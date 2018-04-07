# import sys
# if sys.version_info[0] < 3:
#     import Tkinter as tk
#     from Tkinter import filedialog as fd
# else:
#     import tkinter as tk
#     from tkinter import filedialog as fd

import tkinter as tk
from tkinter import filedialog as fd

import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import gprpy as gp

from scipy import signal
import numpy as np

class GPRPyApp:

    def __init__(self,master):
        # Prepare the frame
        frame = tk.Frame(master)        
        frame.pack()
        self.window = master

        
        # Initialize the gprpy
        proj = gp.gprpy2d()

        # Show splash screen
        fig=Figure(figsize=(8,5))
        a=fig.add_subplot(111)
        splash=signal.ricker(100,4)
        a.plot(splash)
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().pack(side="bottom")
        canvas.draw()
 
        # Initialize plotting variables
        #plotvars = {}
        #plotvars["maxyval"] = 10000000
        
        # Load data
        self.LoadButton = tk.Button(
            frame, text="Import Data", fg="black",
            command=lambda : self.loadData(proj))            
        self.LoadButton.pack(side="left")

        # Refreshing plot
        self.plotButton = tk.Button(
            frame, text="Refresh Plot",
            command=lambda : self.plotTWTTData(proj,fig=fig,a=a,canvas=canvas,
                                               maxyval=float(myv.get()),
                                               contrast=float(contr.get())))
        self.plotButton.pack(side="left")

        # Print Figure
        #self.printButton = tk.Button(
        #    frame, text="Print Figure", fg="black",
        #    command=lambda : self.printFigure(fig=fig,a=a,canvas=canvas))
        
        # Mode switch
        #modeswitch = tk.Listbox(master)
        #modeswitch.pack()
        #modeswitch.insert(tk.END, "plotting mode")
        #for item in ["TWTT", "depth", "WARR", "CMP"]:
        #    modeswitch.insert(tk.END,item)
             
        # y limit
        myvtext = tk.StringVar()
        myvtext.set("Maximum y value")
        myvlabel = tk.Label(master, textvariable=myvtext,height=4)
        myvlabel.pack(side="left")
        myv = tk.StringVar()
        maxybox = tk.Entry(master, textvariable=myv)
        maxybox.pack(side="left")
        myv.set("1000")

        # Contrast
        contrtext = tk.StringVar()
        contrtext.set("Contrast")
        contrlabel = tk.Label(master, textvariable=contrtext,height=4)
        contrlabel.pack(side="left")
        contr = tk.StringVar()
        contrbox = tk.Entry(master, textvariable=contr)
        contrbox.pack(side="left")
        contr.set("1")

        
    def loadData(self,proj):
        filename = fd.askopenfilename( filetypes= (("Sensors and Software", "*.DT1"),
                                                   ("GSSI", ".DZT"),
                                                   ("GPRPy", ".gpr")) )
        #proj = gp.gprpy2d(filename)
        proj.importdata(filename=filename)
        print("Loaded " + filename)


    def plotTWTTData(self,proj,fig,a,canvas,maxyval,contrast):
        print("plotting, y-max " + str(maxyval))
        color="gray"

        stdcont = np.argmax(abs(proj.data))        
        a.imshow(proj.data,cmap=color,extent=[min(proj.profilePos),
                                              max(proj.profilePos),
                                              max(proj.twtt),
                                              min(proj.twtt)],
                 aspect="auto",
                 vmin=-stdcont/contrast, vmax=stdcont/contrast)
        
        #if maxyval < max(proj.twtt):
        a.set_ylim([0,maxyval])
        a.invert_yaxis()

        a.get_xaxis().set_visible(True)
        a.get_yaxis().set_visible(True)
        a.set_ylabel("two-way travel time [ns]")
        a.set_xlabel("profile position")
        a.xaxis.tick_top()
        a.xaxis.set_label_position('top')
        
        # Set contrast
        #if contrast != 1:
        #im=a.get_images()
        #a.set_clim([-stdcont/contrast, stdcont/contrast])
            
        #canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().pack(side="bottom")
        canvas.draw()
        

root = tk.Tk()

app = GPRPyApp(root)

root.mainloop()



