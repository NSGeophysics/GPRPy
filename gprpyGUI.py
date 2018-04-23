# import sys
# if sys.version_info[0] < 3:
#     import Tkinter as tk
#     from Tkinter import filedialog as fd
# else:
#     import tkinter as tk
#     from tkinter import filedialog as fd

import tkinter as tk
from tkinter import filedialog as fd
from tkinter import simpledialog as sd

import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
#import matplotlib.pyplot as plt

import gprpy as gp

from scipy import signal
import numpy as np

colsp=2

class GPRPyApp:

    def __init__(self,master):
        self.window = master

        master.title("GPRPy")
        
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
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=7,rowspan=15,sticky='nsew')
        canvas.draw()
 


        # Load data
        LoadButton = tk.Button(
            text="Import Data", fg="black",
            command=lambda : [self.loadData(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                maxyval=float(myv.get()),
                                                contrast=float(contr.get()),
                                                color=colvar.get())])
        LoadButton.config(height = 1, width = 10)         
        LoadButton.grid(row=0, column=7, sticky='nsew',columnspan=colsp)


        # Dewow
        DewowButton = tk.Button(
            text="Dewow", fg="black",
            command=lambda : [self.dewow(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                maxyval=float(myv.get()),
                                                contrast=float(contr.get()),
                                                color=colvar.get())])
        DewowButton.config(height = 1, width = 10)         
        DewowButton.grid(row=1, column=7, sticky='nsew',columnspan=colsp)


        
        # TimeZero Adjust
        TZAButton = tk.Button(
            text="Time Zero Adj", fg="black",
            command=lambda : [proj.timeZeroAdjust(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                maxyval=float(myv.get()),
                                                contrast=float(contr.get()),
                                                color=colvar.get())])
        TZAButton.config(height = 1, width = 10)         
        TZAButton.grid(row=2, column=7, sticky='nsew',columnspan=colsp)

        
        

        # Rem mean trace
        remMeanTraceButton = tk.Button(
            text="Rem mean tr", fg="black",
            command=lambda : [self.remMeanTrace(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                maxyval=float(myv.get()),
                                                contrast=float(contr.get()),
                                                color=colvar.get())])
        remMeanTraceButton.config(height = 1, width = 10)         
        remMeanTraceButton.grid(row=3, column=7, sticky='nsew',columnspan=colsp)



        # Gain: row 4
        tpowButton = tk.Button(
            text="tpow", fg="black",
            command=lambda : [self.tpowGain(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                   maxyval=float(myv.get()),
                                                   contrast=float(contr.get()),
                                                   color=colvar.get())])
        tpowButton.config(height=1, width=1)
        tpowButton.grid(row=4, column=7, sticky='nsew')

        agcButton = tk.Button(
            text="AGC",fg="black", command=None)
        agcButton.config(height=1, width=1)
        agcButton.grid(row=4, column=8, sticky='nsew')
        

        # Set Velocity: row 11
        setVelButton = tk.Button(
            text="Set velocity", fg="black",
            command=lambda : [self.setVelocity(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                   maxyval=float(myv.get()),
                                                   contrast=float(contr.get()),
                                                   color=colvar.get())])
        setVelButton.config(height = 1, width = 10)         
        setVelButton.grid(row=11, column=7, sticky='nsew',columnspan=colsp)

        # Topo Correct row 12
        


        
        # Save data
        SaveButton = tk.Button(
            text="Save Data", fg="black",
            command=lambda : self.saveData(proj))
        SaveButton.config(height = 1, width = 10)         
        SaveButton.grid(row=13, column=7, sticky='nsew',columnspan=colsp)

        # Print Figure
        PrintButton = tk.Button(
            text="Print Figure", fg="black",
            command=lambda : self.printProfileFig(proj=proj,fig=fig,maxyval=myv.get(),contrast=contr.get(),color=colvar.get()))
        PrintButton.config(height = 1, width = 10)         
        PrintButton.grid(row=14, column=7, sticky='nsew',columnspan=colsp)

        # Write history
        HistButton = tk.Button(
            text="Write history", fg="black",
            command=lambda : self.writeHistory(proj))
        HistButton.config(height = 1, width = 10)         
        HistButton.grid(row=15, column=7, sticky='nsew',columnspan=colsp)
        


        
        ## Plotting
        
        # Refreshing plot
        plotButton = tk.Button(
            text="Refresh Plot",
            command=lambda : self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                               maxyval=myv.get(),
                                               contrast=contr.get(),
                                               color=colvar.get()))
        plotButton.config(height = 1, width = 10)
        plotButton.grid(row=0, column=6, sticky='nsew')

        # Undo Button
        undoButton = tk.Button(
            text="Undo",
            command=lambda : [proj.undo(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas,
                                                maxyval=float(myv.get()),
                                                contrast=float(contr.get()),
                                                color=colvar.get())])
        undoButton.config(height = 1, width = 10)
        undoButton.grid(row=0, column=0, sticky='nsew')
        
                     
             
        # y limit
        myvtext = tk.StringVar()
        myvtext.set("Max y value")
        myvlabel = tk.Label(master, textvariable=myvtext,height = 1,width = 10)
        myvlabel.grid(row=0, column=1, sticky='nsew')
        myv = tk.DoubleVar()
        maxybox = tk.Entry(master, textvariable=myv)
        maxybox.grid(row=0, column=2, sticky='nsew')
        maxybox.config(width=8)
        myv.set("1000000.0")

        # Contrast
        contrtext = tk.StringVar()
        contrtext.set("Contrast")
        contrlabel = tk.Label(master, textvariable=contrtext,height = 1,width = 8)
        contrlabel.grid(row=0, column=3, sticky='nsew')
        contr = tk.DoubleVar()
        contrbox = tk.Entry(master, textvariable=contr, width=8)
        contrbox.grid(row=0, column=4, sticky='nsew')
        contr.set("1.0")

        # Mode switch for figure color
        colvar=tk.StringVar()
        colvar.set("gray")
        colswitch = tk.OptionMenu(master,colvar,"gray","bwr")
        colswitch.grid(row=0, column=5, sticky='nsew')




    def dewow(self,proj):
        window = sd.askinteger("Input","Dewow window width (number of samples)")
        proj.dewow(window=window)
        

    def remMeanTrace(self,proj):
        ntraces = sd.askinteger("Input","Remove mean over how many traces?")
        proj.remMeanTrace(ntraces=ntraces)


    def tpowGain(self,proj):
        power = sd.askfloat("Input","Power for tpow gain?")
        proj.tpowGain(power=power)
        
        
    def setVelocity(self,proj):
        velocity =  sd.askfloat("Input","Radar wave velocity [m/ns]?")
        proj.setVelocity(velocity)

        
    def loadData(self,proj):
        filename = fd.askopenfilename( filetypes= (("GPRPy (.gpr)", "*.gpr"),
                                                   ("Sensors and Software (.DT1)", "*.DT1"),
                                                   ("GSSI (.DZT)", "*.DZT")))
        proj.importdata(filename=filename)
        print("Loaded " + filename)

        
    def saveData(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".gpr")
        proj.save(filename)
       
        
    def writeHistory(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".py")
        proj.writeHistory(filename)
        print("Wrote history to " + filename)


    def plotProfileData(self,proj,fig,a,canvas,maxyval,contrast,color):
        print("plotting, y-max " + str(maxyval))
        #color="gray"

        a.clear()
        
        stdcont = np.argmax(abs(proj.data))
        if proj.velocity is None:
            a.imshow(proj.data,cmap=color,extent=[min(proj.profilePos),
                                                  max(proj.profilePos),
                                                  max(proj.twtt),
                                                  min(proj.twtt)],
                     aspect="auto",
                     vmin=-stdcont/contrast, vmax=stdcont/contrast)
            a.set_ylim([0,min(maxyval,max(proj.twtt))])
            a.set_ylabel("two-way travel time [ns]")
        else:
            a.imshow(proj.data,cmap=color,extent=[min(proj.profilePos),
                                                  max(proj.profilePos),
                                                  max(proj.depth),
                                                  min(proj.depth)],
                     aspect="auto",
                     vmin=-stdcont/contrast, vmax=stdcont/contrast)
            a.set_ylim([0,min(maxyval,max(proj.depth))])
            a.set_ylabel("depth [m]")

        if proj.topoCorrected:
            a.set_ylabel("elevation [m]")
            
        a.invert_yaxis()
        a.get_xaxis().set_visible(True)
        a.get_yaxis().set_visible(True)                    
        a.set_xlabel("profile position")
        a.xaxis.tick_top()
        a.xaxis.set_label_position('top')
        
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=7, rowspan=15, sticky='nsew')
        canvas.draw()
        

    def printProfileFig(self,proj,fig,maxyval,contrast,color):
        figname = fd.asksaveasfilename(defaultextension=".pdf")        
        fig.savefig(figname, format='pdf')        
        # Put what you did in history        
        histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, timelim=[0,%g])" %(figname,color,contrast,maxyval)
        proj.history.append(histstr)
        print("Saved figure as %s" %(figname+'.pdf'))
        
root = tk.Tk()

app = GPRPyApp(root)

root.mainloop()



