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
import toolbox.splash as splash
import os
import Pmw
import scipy.interpolate as interp



colsp=2
rightcol=7
halfwid=4

class GPRPyCWApp:

    def __init__(self,master):
        self.window = master

        master.title("GPRPy CW")

        self.balloon = Pmw.Balloon()
        fig = Figure(figsize=(9,5))
        ahyp = fig.add_axes([0.04,0.1,0.29,0.88])
        adata = fig.add_axes([0.365,0.1,0.29,0.88])
        alin = fig.add_axes([0.69,0.1,0.29,0.88])
        adata.get_yaxis().set_visible(False)
        alin.get_yaxis().set_visible(False)

        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().grid(row=2,column=0,columnspan=7,rowspan=15,sticky='nsew')
        canvas.draw()

        proj = gp.gprpyCW()



        ## Data processing buttons
        
        # Load data
        LoadButton = tk.Button(
            text="import data", fg="black",
            command=lambda : [self.loadData(proj),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        LoadButton.config(height = 1, width = 2*halfwid)         
        LoadButton.grid(row=0, column=rightcol, sticky='nsew',columnspan=colsp,rowspan=2)


        # Adjust profile length; if trigger wheel is not good
        AdjProfileButton = tk.Button(
            text="adj profile", fg="black",
            command=lambda : [self.adjProfile(proj),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        AdjProfileButton.config(height = 1, width = 2*halfwid)         
        AdjProfileButton.grid(row=2, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(AdjProfileButton,
                          "Adjust the profile length to \n"
                          "known start and end positions.")


        # Set new zero time
        SetZeroTimeButton = tk.Button(
            text="set zero time", fg="black",
            command=lambda : [self.setZeroTime(proj),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        SetZeroTimeButton.config(height = 1, width = 2*halfwid)         
        SetZeroTimeButton.grid(row=3, column=rightcol, sticky='nsew',columnspan=colsp)    
        self.balloon.bind(SetZeroTimeButton,
                          "Set the two-way travel time that \n" 
                          "that corresponds to the surface.")

        # truncate Y
        truncYButton = tk.Button(
            text="truncate Y", fg="black",
            command=lambda : [self.truncateY(proj),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        truncYButton.config(height = 1, width = 2*halfwid)         
        truncYButton.grid(row=4, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(truncYButton,
                          "Remove data points at arrival times\n"
                          "later than the chosen value. If velocity\n"
                          "is given: remove data points at depths greater\n"
                          "than the chosen value")   


        # Dewow
        DewowButton = tk.Button(
            text="dewow", fg="black",
            command=lambda : [self.dewow(proj),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        DewowButton.config(height = 1, width = 2*halfwid)         
        DewowButton.grid(row=5, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(DewowButton,
                          "Trace-wise low-cut filter. Removes\n" 
                          "from each trace a running mean of\n"
                          "chosen window width.") 

        
        # Normalize
        NormalizeButton = tk.Button(
            text="normalize", fg="black",
            command = lambda : [proj.normalize(),
                                self.plotCWData(proj,a=adata,canvas=canvas)])
        NormalizeButton.config(height = 1, width = 2*halfwid)         
        NormalizeButton.grid(row=6, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(NormalizeButton,
                          "Normalizes each trace such that\n" 
                          "they all have equal energy") 



                
        # Write history
        HistButton = tk.Button(
            text="write history", fg="black",
            command=lambda : self.writeHistory(proj))
        HistButton.config(height = 1, width = 2*halfwid)         
        HistButton.grid(row=16, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(HistButton,
                          'Saves the history of the current status as a\n'
                          'python script which can be run directly from\n'
                          'the command prompt. If the current data is\n'  
                          'from a .gpr file, then the python script will\n'
                          'contain all steps going back to the raw data.\n'
                          'The script will not contain visualization settings\n'
                          'such as x-range settings, unless the "print figure"\n'
                          'command was used.')




        
        ## Visualization buttons

        # Undo Button
        undoButton = tk.Button(
            text="undo",
            command=lambda : [proj.undo(),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        undoButton.config(height = 1, width = 2*halfwid)
        undoButton.grid(row=0, column=0, sticky='nsew',rowspan=2)
        self.balloon.bind(undoButton,
                          '"Undoes" the most recent processing step and\n' 
                          'sets the data back to its previous state.\n' 
                          'This also removes the most recent processing\n'
                          'step from the history. Does not revert\n' 
                          'visualization settings such as "set x-range"\n'
                          'etc.')

        # X range
        XrngButton = tk.Button(
            text="set x-range", fg="black",
            command=lambda : [self.setXrng(),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        XrngButton.config(height = 1, width = 2*halfwid)         
        XrngButton.grid(row=0, column=1, sticky='nsew',rowspan=2)
        self.balloon.bind(XrngButton,"Set the x-axis display limits.")
        
        
        # Y range
        YrngButton = tk.Button(
            text="set y-range", fg="black",
            command=lambda : [self.setYrng(),
                              self.plotCWData(proj,a=adata,canvas=canvas)])
        YrngButton.config(height = 1, width = 2*halfwid)         
        YrngButton.grid(row=0, column=2, sticky='nsew',rowspan=2)
        self.balloon.bind(YrngButton,"Set the y-axis display limits.")


        
        # Contrast
        contrtext = tk.StringVar()
        contrtext.set("contrast")
        contrlabel = tk.Label(master, textvariable=contrtext,height = 1,width = 2*halfwid)
        contrlabel.grid(row=0, column=4, sticky='nsew')
        self.balloon.bind(contrlabel,"Set color saturation")
        self.contrast = tk.DoubleVar()
        contrbox = tk.Entry(master, textvariable=self.contrast, width=2*halfwid)
        contrbox.grid(row=1, column=4, sticky='nsew')
        #contr.set("1.0")
        self.contrast.set("1.0")

        
        # Mode switch for figure color
        self.color=tk.StringVar()
        self.color.set("gray")
        colswitch = tk.OptionMenu(master,self.color,"gray","bwr")
        colswitch.grid(row=0, column=5, sticky='nsew',rowspan=2)
        self.balloon.bind(colswitch,
                          "Choose between gray-scale\n"
                          "and red-white-blue (rwb)\n" 
                          "data representation.")


        # Refreshing plot
        plotButton = tk.Button(
            text="refresh plot",
            command=lambda : self.plotCWData(proj,a=adata,canvas=canvas))
        plotButton.config(height = 1, width = 2*halfwid)
        plotButton.grid(row=0, column=6, sticky='nsew',rowspan=2)
        self.balloon.bind(plotButton,
                          "Refreshes the figure after changes\n"
                          "in the visualization settings. Also\n"
                          "removes any plotted hyperbolae.")




        
        

    # functions

    
    def loadData(self,proj):
        filename = fd.askopenfilename( filetypes= (("Sensors and Software (.DT1)", "*.DT1"),
                                                   ("GSSI (.DZT)", "*.DZT")))
        if filename is not '':
            self.getType()
            proj.importdata(filename,self.dtype)
            self.xrng = [np.min(proj.profilePos),np.max(proj.profilePos)]
            self.yrng = [0,np.max(proj.twtt)]
           

    def getType(self):
        typeQuery = tk.Toplevel(self.window)
        typeQuery.title("Data type")
        text = tk.Label(typeQuery,text="Is this a common midpoint (CMP) or\n"
                        "wide aperture reflection and refraction (WARR)\n"
                        "data set?",fg="red")
        text.pack(padx=10,pady=10)
        CMPButton = tk.Button(typeQuery, text="CMP", width=10,
                              command = lambda : [self.setCMP(),
                                                  typeQuery.destroy()])
        CMPButton.pack(side="left")
        WARRButton = tk.Button(typeQuery, text="WARR", width=10,
                               command = lambda : [self.setWARR(),
                                                   typeQuery.destroy()])
        WARRButton.pack(side="right")
        # This forces the program to wait until this query is answered
        self.window.wait_window(typeQuery) 
    def setCMP(self):
        self.dtype = 'CMP'
        print("Data type is CMP")
    def setWARR(self):
        self.dtype = 'WARR'
        print("Data type is WARR")


    def plotCWData(self,proj,a,canvas):
        a.clear()
        stdcont = np.nanmax(np.abs(proj.data)[:])
        a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos),
                                                         max(proj.profilePos),
                                                         max(proj.twtt),
                                                         min(proj.twtt)],
                 aspect="auto",
                 vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
        a.set_ylim(self.yrng)
        a.set_xlim(self.xrng)
        a.invert_yaxis()
        if self.dtype is "WARR":
            a.set_xlabel("antenna separation [m]")
        elif self.dtype is "CMP":
            a.set_xlabel("distance from midpoint [m]")
        # Allow for cursor coordinates being displayed        
        def moved(event):
            if event.xdata is not None and event.ydata is not None:
                canvas.get_tk_widget().itemconfigure(tag, text="(x = %5.5g, y = %5.5g)" % (event.xdata, event.ydata))
                
        canvas.mpl_connect('button_press_event', moved)
        tag = canvas.get_tk_widget().create_text(20, 20, text="", anchor="nw")

        canvas.get_tk_widget().grid(row=2,column=0,columnspan=7, rowspan=15, sticky='nsew')
        canvas.draw()



    def setYrng(self):
        ylow = sd.askfloat("Input","Min Y value")
        if ylow is not None:            
            yhigh = sd.askfloat("Input","Max Y value")
            if yhigh is not None:
                self.prevyrng=self.yrng
                self.yrng=[ylow,yhigh]


    def setXrng(self):
        xlow = sd.askfloat("Input","Min X value")
        if xlow is not None:
            xhigh = sd.askfloat("Input","Max X value")
            if xhigh is not None:
                self.xrng=[xlow,xhigh]
        

    def adjProfile(self,proj):
        minPos = sd.askfloat("Input","Start x coordinate")
        if minPos is not None:
            maxPos = sd.askfloat("Input","End x coordinate")
            if maxPos is not None:
                proj.adjProfile(minPos=minPos,maxPos=maxPos)
                self.xrng=[minPos,maxPos]

    def setZeroTime(self,proj):
        newZeroTime = sd.askfloat("Input","New zero time")
        if newZeroTime is not None:
            proj.setZeroTime(newZeroTime=newZeroTime)
                


    def writeHistory(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".py")
        if filename is not '':
            proj.writeHistory(filename)
            print("Wrote history to " + filename)
            
            
         
root = tk.Tk()

root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

app = GPRPyCWApp(root)

root.mainloop()
