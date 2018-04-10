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
        self.window = master

        
        # Initialize the gprpy
        proj = gp.gprpy2d()

        # Show splash screen
        fig=Figure(figsize=(8,5))
        a=fig.add_subplot(111)
        splash=signal.ricker(50,4)
        a.plot(splash)
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        #canvas.get_tk_widget().pack(side="bottom")
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=6,rowspan=8)
        canvas.draw()
 
        # Initialize plotting variables
        #plotvars = {}
        #plotvars["maxyval"] = 10000000
        

        # Refreshing plot
        plotButton = tk.Button(
            text="Refresh Plot",
            command=lambda : self.plotTWTTData(proj,fig=fig,a=a,canvas=canvas,
                                               maxyval=float(myv.get()),
                                               contrast=float(contr.get())))
        plotButton.config(height = 2, width = 10)
        plotButton.grid(row=0,column=0)

        # Undo Button
        undoButton = tk.Button(
            text="Undo",
            command=lambda : proj.undo())
        undoButton.config(height = 2, width = 10)
        undoButton.grid(row=0,column=1)
        
        # Load data
        LoadButton = tk.Button(
            text="Import Data", fg="black",
            command=lambda : self.loadData(proj))
        LoadButton.config(height = 2, width = 10)         
        LoadButton.grid(row=0,column=6)

        # TimeZero Adjust
        TZAButton = tk.Button(
            text="Time Zero Adj", fg="black",
            command=lambda : proj.timeZeroAdjust())
        TZAButton.config(height = 2, width = 10)         
        TZAButton.grid(row=1,column=6)
        
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
        myvtext.set("Max y value")
        myvlabel = tk.Label(master, textvariable=myvtext,height = 2,width = 8)
        myvlabel.grid(row=0,column=2)

        myv = tk.StringVar()
        maxybox = tk.Entry(master, textvariable=myv)
        maxybox.grid(row=0,column=3)
        maxybox.config(width=8)
        myv.set("1000000")

        # Contrast
        contrtext = tk.StringVar()
        contrtext.set("Contrast")
        contrlabel = tk.Label(master, textvariable=contrtext,height = 2,width = 8)
        contrlabel.grid(row=0,column=4)
        contr = tk.StringVar()
        contrbox = tk.Entry(master, textvariable=contr,width=8)
        contrbox.grid(row=0,column=5)
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
        a.set_ylim([0,min(maxyval,max(proj.twtt))])
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
        #canvas.get_tk_widget().pack(side="bottom")
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=5)
        #canvas.get_tk_widget().place(relx=0.,rely=0.2,anchor="c")
        canvas.draw()
        

root = tk.Tk()

app = GPRPyApp(root)

root.mainloop()



