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


colsp=2
rightcol=7
halfwid=6


class GPRPyApp:

    def __init__(self,master):
        self.window = master

        master.title("GPRPy")
        
        # Variables specific to GUI
        self.balloon = Pmw.Balloon()
        self.picking = False
        
        # Initialize the gprpy
        proj = gp.gprpy2d()

        # Show splash screen
        fig=Figure(figsize=(8,5))
        a=fig.add_subplot(111)
        dir_path = os.path.dirname(os.path.realpath(__file__))
        splash.showSplash(a,dir_path)
        
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().grid(row=2,column=0,columnspan=7,rowspan=15,sticky='nsew')
        canvas.draw() 



        ## Visualization Buttons              

        # Undo Button
        undoButton = tk.Button(
            text="undo",
            command=lambda : [self.resetYrng(proj),
                              proj.undo(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
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
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        XrngButton.config(height = 1, width = 2*halfwid)         
        XrngButton.grid(row=0, column=1, sticky='nsew',rowspan=2)
        self.balloon.bind(XrngButton,"Set the x-axis display limits.")
        

        # Y range
        YrngButton = tk.Button(
            text="set y-range", fg="black",
            command=lambda : [self.setYrng(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        YrngButton.config(height = 1, width = 2*halfwid)         
        YrngButton.grid(row=0, column=2, sticky='nsew',rowspan=2)
        self.balloon.bind(YrngButton,"Set the y-axis display limits.")

        # Aspect
        AspButton = tk.Button(
            text="aspect ratio", fg="black",
            command=lambda : [self.setAspect(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])                              
        AspButton.config(height = 1, width = 2*halfwid)         
        AspButton.grid(row=0, column=3, sticky='nsew',rowspan=2)
        self.balloon.bind(AspButton, "Set the aspect ratio between x- and y-axis.")
        

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
            command=lambda : self.plotProfileData(proj,fig=fig,a=a,canvas=canvas))
        plotButton.config(height = 1, width = 2*halfwid)
        plotButton.grid(row=0, column=6, sticky='nsew',rowspan=2)
        self.balloon.bind(plotButton,
                          "Refreshes the figure after changes\n"
                          "in the visualization settings. Also\n"
                          "removes any plotted hyperbolae.")




        ## Methods buttons
        
        
        # Load data
        LoadButton = tk.Button(
            text="import data", fg="black",
            command=lambda : [self.loadData(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        LoadButton.config(height = 1, width = 2*halfwid)         
        LoadButton.grid(row=0, column=rightcol, sticky='nsew',columnspan=colsp,rowspan=2)
        self.balloon.bind(LoadButton,"Load .gpr, .DT1, or .DZT data.")

        # Adjust profile length; if trigger wheel is not good
        AdjProfileButton = tk.Button(
            text="adj profile", fg="black",
            command=lambda : [self.adjProfile(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        AdjProfileButton.config(height = 1, width = 2*halfwid)         
        AdjProfileButton.grid(row=2, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(AdjProfileButton,
                          "Adjust the profile length to \n"
                          "known start and end positions.")

        # Set new zero time
        SetZeroTimeButton = tk.Button(
            text="set zero time", fg="black",
            command=lambda : [self.setZeroTime(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        SetZeroTimeButton.config(height = 1, width = 2*halfwid)         
        SetZeroTimeButton.grid(row=3, column=rightcol, sticky='nsew',columnspan=colsp)    
        self.balloon.bind(SetZeroTimeButton,
                          "Set the two-way travel time that \n" 
                          "that corresponds to the surface.")
                

        # truncate Y
        truncYButton = tk.Button(
            text="truncate Y", fg="black",
            command=lambda : [self.truncateY(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
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
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        DewowButton.config(height = 1, width = 2*halfwid)         
        DewowButton.grid(row=5, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(DewowButton,
                          "Trace-wise low-cut filter. Removes\n" 
                          "from each trace a running mean of\n"
                          "chosen window width.")                          

        
        # TimeZero Adjust
        TZAButton = tk.Button(
            text="time zero adj", fg="black",
            command=lambda : [proj.timeZeroAdjust(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        TZAButton.config(height = 1, width = 2*halfwid)         
        TZAButton.grid(row=6, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(TZAButton,
                         'Automatically shifts each trace up or down\n'
                         'such that the maximum aplitudes of the individual\n'
                         'traces align. Can lead to problems when the maxima\n' 
                         'are not in the air waves. If the results are bad,\n' 
                         'use the "undo" button.')
        
        

        # Rem mean trace
        remMeanTraceButton = tk.Button(
            text="rem mean tr", fg="black",
            command=lambda : [self.remMeanTrace(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        remMeanTraceButton.config(height = 1, width = 2*halfwid)         
        remMeanTraceButton.grid(row=7, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(remMeanTraceButton,
                          "Removes from each traces the average\n" 
                          "of its surrounding traces. This can be\n"
                          "useful to remove air waves, or\n" 
                          "horizontal features.")


        # Gain: row 7
        tpowButton = tk.Button(
            text="tpow", fg="black",
            command=lambda : [self.tpowGain(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        tpowButton.config(height=1, width=halfwid)
        tpowButton.grid(row=8, column=rightcol, sticky='nsew')
        self.balloon.bind(tpowButton,
                          "t-power gain. Increases the power of the\n"
                          "signal by a factor of (two-way travel time)^p,\n"
                          "where the user provides p. This gain is often\n" 
                          "less aggressive than agc.")

        
        agcButton = tk.Button(
            text="agc",fg="black",
            command=lambda : [self.agcGain(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        agcButton.config(height=1, width=halfwid)
        agcButton.grid(row=8, column=rightcol+1, sticky='nsew')
        self.balloon.bind(agcButton,
                          "Automatic gain controll. Normalizes the power\n"
                          "of the signal per given sample window along\n" 
                          "each trace.")

        # show hyperbola
        hypButton = tk.Button(
            text="show hyperb", fg="black",
            command=lambda : [self.showHyp(proj,a), canvas.draw()])
        hypButton.config(height = 1, width = 2*halfwid)
        hypButton.grid(row=9, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(hypButton,
                          "Draws a hyperbola depending on profile position,\n"
                          "two-way travel time, and estimated velocity. This\n" 
                          "can be used to find the subsurface velocity when\n"
                          "a hyperbola is visible in the data.\n"
                          "The plotted hyperbola will disappear when the image\n" 
                          "is refreshed.")

        

        # Set Velocity: row 8
        setVelButton = tk.Button(
            text="set velocity", fg="black",
            command=lambda : [self.setVelocity(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        setVelButton.config(height = 1, width = 2*halfwid)         
        setVelButton.grid(row=10, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(setVelButton,
                          "Set the known subsurface radar velocity. This will\n" 
                          "turn the y-axis from two-way travel time to depth.\n"
                          "This step is necessary for topographic correction.")

        # Topo Correct row 9
        topoCorrectButton = tk.Button(
            text="topo correct", fg="black",
            command=lambda : [self.topoCorrect(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        topoCorrectButton.config(height = 1, width = 2*halfwid)
        topoCorrectButton.grid(row=11, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(topoCorrectButton,
                          "Reads a comma- or tab-separated file containing\n" 
                          "either 3 columns (easting, northing, elevation)\n" 
                          "or two columns (profile position, elevation).\n" 
                          "All coordinates in meters.")                                                      


       
        startPickButton = tk.Button(
            text="start pick", fg="black",
            command=lambda : self.startPicking(proj,fig=fig,a=a,canvas=canvas))        
        startPickButton.config(height = 1, width = halfwid)
        startPickButton.grid(row=12, column=rightcol, sticky='nsew',columnspan=1)


        stopPickButton = tk.Button(
            text="stop pick", fg="black",
            command=lambda : [self.stopPicking(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        
        stopPickButton.config(height = 1, width = halfwid)
        stopPickButton.grid(row=12, column=rightcol+1, sticky='nsew',columnspan=1)
        
        
        # Save data
        SaveButton = tk.Button(
            text="save data", fg="black",
            command=lambda : self.saveData(proj))
        SaveButton.config(height = 1, width = 2*halfwid)         
        SaveButton.grid(row=13, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(SaveButton,
                          'saves the processed data including its history in a\n'
                          '.gpr file. The resulting file will contain absolute\n'
                          'path names of the used data and topography files.\n'
                          'Visualization settings such as "set x-range" or\n'
                          '"contrast" will not be saved.')

        
        
        # Print Figure
        PrintButton = tk.Button(
            text="print figure", fg="black",
            command=lambda : self.printProfileFig(proj=proj,fig=fig))
        PrintButton.config(height = 1, width = 2*halfwid)         
        PrintButton.grid(row=14, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(PrintButton,
                          "Saves the current visible figure in a pdf with \n"
                          "chosen resolution. If there is a hyperbola on\n" 
                          "the current figure, then the hyperbola will also\n"
                          "appear on the printed figure.")


        # Export to VTK
        VTKButton = tk.Button(
            text="export to VTK", fg="black",
            command = lambda : self.exportVTK(proj))
        VTKButton.config(height = 1, width = 2*halfwid)
        VTKButton.grid(row=15, column=rightcol, sticky='nsew',columnspan=colsp)
        self.balloon.bind(VTKButton,
                          "Exports the processed figure to a\n"
                          "VTK format, that can be read by\n" 
                          "Paraview or similar 3D programs")
        


        
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


        
      
    def setYrng(self):
        ylow = sd.askfloat("Input","Min Y value")
        yhigh = sd.askfloat("Input","Max Y value")
        self.prevyrng=self.yrng
        self.yrng=[ylow,yhigh]
        

    def resetYrng(self,proj):
        # Only needed in undo, and only if what you want to
        # undo changed the y axis
        if ("setVelocity" in proj.history[-1]) or ("topoCorrect" in proj.history[-1]): 
            self.yrng=self.prevyrng


    def setAspect(self):
        self.asp = sd.askfloat("Input","Plotting aspect ratio")
        

    def setXrng(self):
        xlow = sd.askfloat("Input","Min X value")
        xhigh = sd.askfloat("Input","Max X value")
        self.xrng=[xlow,xhigh]
        

    def adjProfile(self,proj):
        minPos = sd.askfloat("Input","Start x coordinate")
        maxPos = sd.askfloat("Input","End x coordinate")
        proj.adjProfile(minPos=minPos,maxPos=maxPos)
        self.xrng=[minPos,maxPos]

    def setZeroTime(self,proj):
        newZeroTime = sd.askfloat("Input","New zero time")
        proj.setZeroTime(newZeroTime=newZeroTime)
        
        
    def dewow(self,proj):
        window = sd.askinteger("Input","Dewow window width (number of samples)")
        proj.dewow(window=window)
        

    def remMeanTrace(self,proj):
        ntraces = sd.askinteger("Input","Remove mean over how many traces?")
        proj.remMeanTrace(ntraces=ntraces)


    def tpowGain(self,proj):
        power = sd.askfloat("Input","Power for tpow gain?")
        proj.tpowGain(power=power)
        

    def agcGain(self,proj):
        window = sd.askinteger("Input","Window length for AGC?")
        proj.agcGain(window=window)

    def truncateY(self,proj):
        maxY = sd.askfloat("Input","Truncate at what y value\n" 
                           "(two-way travel time or depth)")
        proj.truncateY(maxY)
        
        
    def setVelocity(self,proj):
        velocity =  sd.askfloat("Input","Radar wave velocity [m/ns]?")
        proj.setVelocity(velocity)
        self.prevyrng=self.yrng
        self.yrng=[0,np.max(proj.depth)]

    def topoCorrect(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Topo Correct Error","You have to set the velocity first")
            return
        topofile = fd.askopenfilename()
        commasep = mesbox.askyesno("Question","Is this a comma-separated file (Yes)\nor tab-separated (No)")
        if commasep:
            delimiter = ','
        else:
            delimiter = '\t'            
        proj.topoCorrect(topofile,delimiter)
        self.prevyrng=self.yrng
        self.yrng=[proj.maxTopo-np.max(proj.depth),proj.maxTopo]



   

    def startPicking(self,proj,fig,a,canvas):
        self.picking = True
        self.picked = np.asmatrix(np.empty((0,2)))
        print("Picking mode on")
        def addPoint(event):
            self.picked = np.append(self.picked,np.asmatrix([event.xdata,event.ydata]),axis=0)
            self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)
            print('freshly picked')
        canvas.mpl_connect('button_press_event', addPoint)



            

    def stopPicking(self):
        self.picking = False
        print("Picking mode off")
        filename = fd.asksaveasfilename(defaultextension=".txt")
        np.savetxt(filename,self.picked,delimiter='\t')
        print('saved picked file as "%s"' %(filename))
        
        
    def loadData(self,proj):
        filename = fd.askopenfilename( filetypes= (("GPRPy (.gpr)", "*.gpr"),
                                                   ("Sensors and Software (.DT1)", "*.DT1"),
                                                   ("GSSI (.DZT)", "*.DZT")))
        proj.importdata(filename=filename)
        self.xrng = [np.min(proj.profilePos),np.max(proj.profilePos)]
        
        if proj.depth is None:
            self.yrng = [0,np.max(proj.twtt)]
        else:
            if proj.maxTopo is None:
                self.yrng = [0,np.max(proj.depth)]
            else:
                self.yrng = [proj.maxTopo-np.max(proj.depth), proj.maxTopo]

        self.asp=None

        # Just in case someone presses undo before changing yrange        
        self.prevyrng=self.yrng    
        print("Loaded " + filename)

        
    def saveData(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".gpr")
        proj.save(filename)


    def exportVTK(self,proj):                    
        outfile = fd.asksaveasfilename()
        gpyes = mesbox.askyesno("Question",
                                "Do you have an x,y,z coordinate file for this profile?")
        delimiter=','
        if gpyes:
            gpsfile = fd.askopenfilename()
            commasep = mesbox.askyesno("Question","Is this a comma-separated file (Yes)\nor tab-separated (No)")
            if commasep:
                delimiter = ','
            else:
                delimiter = '\t' 
        else: 
            gpsfile = None
        
        thickness = sd.askfloat("Input","Profile thickness [m]")
        
        if self.asp is None:
            aspect = 1.0
        else:
            aspect = self.asp
            
        proj.exportVTK(outfile,gpsfile=gpsfile,thickness=thickness,delimiter=delimiter,aspect=aspect)
        print('... done with exporting to VTK.')
                
    def writeHistory(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".py")
        proj.writeHistory(filename)
        print("Wrote history to " + filename)


    def plotProfileData(self,proj,fig,a,canvas):
        a.clear()        
        stdcont = np.nanmax(np.abs(proj.data)[:])        
        if proj.velocity is None:
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos),
                                                  max(proj.profilePos),
                                                  max(proj.twtt),
                                                  min(proj.twtt)],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)
            a.set_ylabel("two-way travel time [ns]")
            a.invert_yaxis()
        elif proj.maxTopo is None:
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos),
                                                  max(proj.profilePos),
                                                  max(proj.depth),
                                                  min(proj.depth)],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylabel("depth [m]")
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)
            a.invert_yaxis()
        else:
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos),
                                                  max(proj.profilePos),
                                                  proj.maxTopo-max(proj.depth),
                                                  proj.maxTopo-min(proj.depth)],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylabel("elevation [m]")
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)

        a.get_xaxis().set_visible(True)
        a.get_yaxis().set_visible(True)                    
        a.set_xlabel("profile position [m]")
        a.xaxis.tick_top()
        a.xaxis.set_label_position('top')
        if self.asp is not None:
            a.set_aspect(self.asp)

        # In case you are picking
        if self.picking:
            a.plot(self.picked[:,0],self.picked[:,1],'-x',color='yellow',linewidth=3) 
            a.plot(self.picked[:,0],self.picked[:,1],'-x',color=[0.3,0.5,0.8],linewidth=1)                               

        # Allow for cursor coordinates being displayed        
        def moved(event):
            if event.xdata is not None and event.ydata is not None:
                canvas.get_tk_widget().itemconfigure(tag, text="(x = %5.5g, y = %5.5g)" % (event.xdata, event.ydata))
                
        canvas.mpl_connect('button_press_event', moved)
        tag = canvas.get_tk_widget().create_text(20, 20, text="", anchor="nw")

        canvas.get_tk_widget().grid(row=2,column=0,columnspan=7, rowspan=15, sticky='nsew')
        canvas.draw()


        

    # Show hyperbola
    def showHyp(self,proj,a):
        x0 = sd.askfloat("Input","Hyperbola center on profile [m]")
        t0 = sd.askfloat("Input","Hyperbola apex location (two-way travel time [ns])")
        v  = sd.askfloat("Input","Estimated velocity [m/ns]")
        # t0 is two-way travel time
        y=proj.profilePos-x0
        d=v*t0/2.0
        k=np.sqrt(d**2 + np.power(y,2))
        t2=2*k/v
        a.plot(proj.profilePos,t2,'--c',linewidth=3)
        

    def printProfileFig(self,proj,fig):
        figname = fd.asksaveasfilename(defaultextension=".pdf")
        dpi = sd.askinteger("Input","Resolution in dots per inch? (Recommended: 600)")
        fig.savefig(figname, format='pdf', dpi=dpi)        
        # Put what you did in history        
        histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], asp=%g, dpi=%d)" %(figname,self.color.get(),self.contrast.get(),self.yrng[0],self.yrng[1],self.xrng[0],self.xrng[1],self.asp,dpi)
        proj.history.append(histstr)
        print("Saved figure as %s" %(figname+'.pdf'))
    

        
root = tk.Tk()

app = GPRPyApp(root)

root.mainloop()




