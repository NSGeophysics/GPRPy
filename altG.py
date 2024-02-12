import sys
from tkinter import ttk
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import simpledialog as sd
from tkinter import messagebox as mesbox
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import gprpy.gprpy as gp
import numpy as np
import gprpy.toolbox.splash as splash
import os
import Pmw
import scipy.interpolate as interp
from tkinter.ttk import *

# Importing Collapsible Pane class that we have
# created in separate file
from collapsiblepane import CollapsiblePane as cp
 
# Making root window or parent window


# Creating Object of Collapsible Pane Container
# If we do not pass these strings in
# parameter the default strings will appear
# on button that were, expand >>, collapse <<
PAD = 4
WIDTH = 10
HEIGHT = 1
HEADING = ('TkDefaultFont', 13,'bold')
BTN_GO_W = 3

class GPRPyApp:
    def __init__(self, master):
        self.window = master

        normscrwidt=1280 #1024
        normscrhigt=720 #768
        self.widfac=normscrwidt/normscrhigt
        self.highfac=1
        #fontfac=(normscrwidt/normscrhigt)/(scrwidt/scrhigt)
        fontfac=1
         # Set some default choices
        self.hypx = 0
        self.hypt = 0
        self.hypv = 0.1

        # Variables specific to GUI
        self.balloon = Pmw.Balloon()
        self.picking = False       
        self.delimiter = None
        self.grid = False

        master.geometry('1024x768')
        master.title("GPR - Py")
        master.rowconfigure(0, minsize=766, weight=1)
        master.columnconfigure(1, minsize=766, weight=1)
       

        proj = gp.gprpyProfile()

        btn_frm = tk.Frame(master, relief= tk.RAISED, bd = 2)
        btn_frm.grid(row = 0, column= 0, sticky= 'ns', rowspan = 2)
        master.rowconfigure(1, weight=10) 

        fig=Figure(figsize=(self.widfac,self.highfac))
        a=fig.add_subplot(111)
        mpl.rcParams.update({'font.size': mpl.rcParams['font.size']*self.widfac})
        a.tick_params(direction='out',length=6*self.widfac,width=self.highfac)
        
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().grid(row=0,column=1,columnspan= 9,rowspan= 22,sticky='nsew')
        #canvas.get_tk_widget().grid(row=2,column=0,columnspan= figcolsp,rowspan= figrowsp,sticky='nsew')
        canvas.draw() 
        
####################################################################################### 
        #File Operation Pane
        #Contains all buttons and labels relation to File Operations
        FM_cpane = cp(btn_frm, 'File Controls -', 'File Controls +')
        FM_cpane.grid(row = 0, column = 0, sticky = 'ew')

        LoadButton = tk.Button(FM_cpane.frame,
            text="Import Data", fg="black",
            command=lambda : [self.loadData(proj), 
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        LoadButton.config(height = HEIGHT, width = WIDTH)         
        LoadButton.grid(row=0, column=0, sticky='nsew',pady = PAD)
        self.balloon.bind(LoadButton,"Load .gpr, .DT1, or .DZT data.")

        undoButton = tk.Button(FM_cpane.frame,
            text="Undo",
            command=lambda : [self.resetYrng(proj),
                              self.undo(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        undoButton.config(height = HEIGHT, width = WIDTH)
        undoButton.grid(row=1, column=0, sticky='ew',pady= PAD )
        self.balloon.bind(undoButton,
                          '"Undoes" the most recent processing step and\n' 
                          'sets the data back to its previous state.\n' 
                          'This also removes the most recent processing\n'
                          'step from the history. Does not revert\n' 
                          'visualization settings such as "set x-range"\n'
                          'etc.')

        

        SaveButton = tk.Button(FM_cpane.frame, 
            text="Save Data", fg="black",
            command=lambda : self.saveData(proj))
        SaveButton.config(height = HEIGHT, width = WIDTH)         
        SaveButton.grid(row= 4, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(SaveButton,
                          'saves the processed data including its history in a\n'
                          '.gpr file. The resulting file will contain absolute\n'
                          'path names of the used data and topography files.\n'
                          'Visualization settings such as "set x-range" or\n'
                          '"contrast" will not be saved.')

        PrintButton = tk.Button(FM_cpane.frame, 
            text="Print Figure", fg="black",
            command=lambda : self.printProfileFig(proj=proj,fig=fig))
        PrintButton.config(height = HEIGHT, width = WIDTH)         
        PrintButton.grid(row=5, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(PrintButton,
                          "Saves the current visible figure in a pdf with \n"
                          "chosen resolution. If there is a hyperbola on\n" 
                          "the current figure, then the hyperbola will also\n"
                          "appear on the printed figure.")


        # Export to VTK
        VTKButton = tk.Button(FM_cpane.frame, 
            text="Export to VTK", fg="black",
            command = lambda : self.exportVTK(proj))
        VTKButton.config(height = HEIGHT, width = WIDTH)
        VTKButton.grid(row=6, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(VTKButton,
                          "Exports the processed figure to a\n"
                          "VTK format, that can be read by\n" 
                          "Paraview or similar 3D programs.")
        


        
        # Write script
        HistButton = tk.Button(FM_cpane.frame, 
            text="Write Script", fg="black",
            command=lambda : self.writeHistory(proj))
        HistButton.config(height = HEIGHT, width = WIDTH)         
        HistButton.grid(row=6, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(HistButton,
                          'Writes a python script to reproduce the \n'
                          'current status.\n'
                          '\n'
                          'If the current data is from a .gpr file, \n'  
                          'then the python script will contain all \n'
                          'steps going back to the raw data. \n'
                          '\n'
                          'The script will not contain visualization \n'
                          'settings such as x-range settings, unless \n'
                          'the "print figure" command was used. ')

####################################################################################### 
        #View Control Pane and all buttons within it
        VC_cpane = cp(btn_frm, 'View Controls -', 'View Controls +')
        VC_cpane.grid(row = 1, column = 0, sticky = 'ew',) 


        # Full view
        FullButton = tk.Button(VC_cpane.frame,
            text="Full View", fg="black",
            command=lambda : [self.setFullView(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        FullButton.config(height = HEIGHT, width = WIDTH)         
        FullButton.grid(row=0, column=0, sticky='ew', padx = PAD, pady = PAD)
        self.balloon.bind(FullButton,"Resets x- and y-axis limits to full data.")

        
        # Grid button
        GridButton = tk.Button(VC_cpane.frame,
            text="Grid", fg="black",
            command=lambda : [self.toggleGrid(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        GridButton.config(height = HEIGHT, width = WIDTH)         
        GridButton.grid(row = 0, column=1, sticky='ew', padx = PAD, pady = PAD)
        self.balloon.bind(GridButton,"Toggles grid on/off.")

        startPickButton = tk.Button(VC_cpane.frame,
            text="start pick", fg="black",
            command=lambda : self.startPicking(proj,fig=fig,a=a,canvas=canvas))        
        startPickButton.config(height = HEIGHT, width = WIDTH) 
        startPickButton.grid(row=1, column=0, sticky='ew',padx = PAD, pady = PAD)
        self.balloon.bind(startPickButton,
                          "Start collecting location information\n" 
                          "by clicking on the profile.")  
        

        stopPickButton = tk.Button(VC_cpane.frame,
            text="stop pick", fg="black",
            command=lambda : [self.stopPicking(proj,canvas),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        stopPickButton.config(height = HEIGHT, width = WIDTH) 
        stopPickButton.grid(row=1, column=1, sticky='ew',padx = PAD, pady = PAD)
        self.balloon.bind(stopPickButton,
                          "Stop collecting location information\n"
                          "and save the locations you collected\n"
                          "in a text file.")
        #Set X-Range
        lbl_xrng = tk.Label(VC_cpane.frame, text= "Set X - Range", font = HEADING)
        lbl_xrng.grid(row = 2, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        
        lbl_minx = tk.Label(VC_cpane.frame, text= "Min:")
        lbl_minx.grid(row = 3, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
    
        tb_minx = tk.Text(VC_cpane.frame, height = 1, width= 3)
        tb_minx.grid(row = 3, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_maxx = tk.Label(VC_cpane.frame, text= "Max:")
        lbl_maxx.grid(row = 3, column = 2, sticky= 'ew', padx= PAD, pady = PAD)

        tb_maxx = tk.Text(VC_cpane.frame, height = 1, width= 3)
        tb_maxx.grid(row = 3, column = 3, sticky= 'ew', padx= PAD, pady = PAD)

        XrngButton = tk.Button(VC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.setXrng(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        XrngButton.config(height = HEIGHT, width = 3)      
        XrngButton.grid(row=3, column=4, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(XrngButton,"Set the x-axis display limits.")
        

        # Y range
        lbl_yrng = tk.Label(VC_cpane.frame, text= "Set Y - Range", font = HEADING)
        lbl_yrng.grid(row = 4, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        
        lbl_miny = tk.Label(VC_cpane.frame, text= "Min:")
        lbl_miny.grid(row = 5, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
    
        tb_miny = tk.Text(VC_cpane.frame, height = 1, width= 3)
        tb_miny.grid(row = 5, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        lbl_maxy = tk.Label(VC_cpane.frame, text= "Max:")
        lbl_maxy.grid(row = 5, column = 2, sticky= 'ew', padx= PAD, pady = PAD)

        tb_maxy = tk.Text(VC_cpane.frame, height = 1, width= 3)
        tb_maxy.grid(row = 5, column = 3, sticky= 'ew', padx= PAD, pady = PAD)

        YrngButton = tk.Button(VC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.setYrng(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        YrngButton.config(height = HEIGHT, width = 3)         
        YrngButton.grid(row=5, column= 4, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(YrngButton,"Set the y-axis display limits.")

        # Contrast
        
        contrlabel = tk.Label(VC_cpane.frame, text  = "Contrast",height = 1,width = WIDTH, font = HEADING)
        contrlabel.grid(row=6, column=0, sticky='nsew')
        self.balloon.bind(contrlabel,"Set color saturation")
        self.contrast = tk.DoubleVar(VC_cpane.frame)
        contrbox = tk.Entry(VC_cpane.frame, textvariable=self.contrast, width=WIDTH)
        contrbox.grid(row=6, column=1, sticky='nsew')
        #contr.set("1.0")
        self.contrast.set("1.0")

        
        # Mode switch for figure color
        self.color=tk.StringVar(VC_cpane.frame)
        self.color.set("gray")
        colswitch = tk.OptionMenu(VC_cpane.frame,self.color,"gray","bwr")
        colswitch.grid(row=6, column=2, sticky='nsew')
        self.balloon.bind(colswitch,
                          "Choose between gray-scale\n"
                          "and blue-white-red (bwr)\n" 
                          "data representation.")
        
    
        # Aspect
        lbl_aspr = tk.Label(VC_cpane.frame, text= "Aspect Ratio", font = HEADING)
        lbl_aspr.grid(row = 7, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_aspr = tk.Text(VC_cpane.frame, height = 1, width= 3)
        tb_aspr.grid(row = 7, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        AspButton = tk.Button(VC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.setAspect(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])                              
        AspButton.config(height = HEIGHT, width = 3)      
        AspButton.grid(row=7, column=3, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(lbl_aspr, "Set the aspect ratio between x- and y-axis.")
    

        #View Control Pane and all buttons within it
        Velo_cpane = cp(btn_frm, 'Velocity Controls -', 'Velocity Controls +')
        Velo_cpane.grid(row = 2, column = 0, sticky = 'ew',) 

        # Set Velocity
        lbl_velo = tk.Label(Velo_cpane.frame, text= "Set Velocity", font = HEADING)
        lbl_velo.grid(row = 0, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_velo = tk.Text(Velo_cpane.frame, height = 1, width= WIDTH)
        tb_velo.grid(row = 0, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        setVelButton = tk.Button(Velo_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.setVelocity(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        setVelButton.config(height = 1, width = 3)         
        setVelButton.grid(row=0, column=3, sticky='nsew')
        self.balloon.bind(setVelButton,
                          "Set the known subsurface radar velocity. This\n" 
                          "turns the y-axis from travel time to depth.\n"
                          "This step is necessary for topographic correction.")


        
        # Correct for antenna separation
        antennaSepButton = tk.Button(Velo_cpane.frame,
            text="Antenna Sep Correct", fg="black",
            command=lambda : [self.antennaSep(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        antennaSepButton.config(height = HEIGHT, width = WIDTH+10)         
        antennaSepButton.grid(row=1, column=0, sticky='ew', padx= PAD, pady = PAD)
        self.balloon.bind(antennaSepButton,                        
                          "If the antenna offset is provided, this corrects for\n"
                          "distortion of arrival times near the surface due to\n"
                          "the separation of transmitter and receiver antenna.\n"
                          "You must have picked the first break of the airwave\n"
                          "for this to function properly and the velocity must be set.")
        

        # Migration Button
        lbl_mig = tk.Label(Velo_cpane.frame, text= "Fk Migration", font = HEADING)
        lbl_mig.grid(row = 2, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_mig = tk.Text(Velo_cpane.frame, height = 1, width= WIDTH )
        tb_mig.grid(row = 2, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        migButton = tk.Button(Velo_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.fkMigration(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        migButton.config(height = 1, width = 3)         
        migButton.grid(row=2, column=3, sticky='nsew')
        self.balloon.bind(migButton,
                          "Stolt's fk migration using a code originally written\n"
                          "in Matlab for the CREWES software package.\n" 
                          "Translated into Python 2 by Nat Wilson.")     

        
        
        # Topo Correct
        lbl_topo = tk.Label(Velo_cpane.frame, text= "Topo Correction", font = HEADING)
        lbl_topo.grid(row = 3, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_topo = tk.Text(Velo_cpane.frame, height = 1, width= WIDTH)
        tb_topo.grid(row = 3, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        topoCorrectButton = tk.Button(Velo_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.topoCorrect(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        topoCorrectButton.config(height = 1, width = 3)
        topoCorrectButton.grid(row=3, column=3, sticky='nsew')
        self.balloon.bind(topoCorrectButton,
                          "Reads a comma- or tab-separated file containing\n" 
                          "either 3 columns (easting, northing, elevation)\n" 
                          "or two columns (profile position, elevation).\n" 
                          "All coordinates in meters.")       

####################################################################################### 
        #Profile Controls
        ProfC_cpane = cp(btn_frm, 'Profile Controls -', 'Profile Controls +')
        ProfC_cpane.grid(row = 3, column = 0, sticky = 'ew',) 

        # Adjust profile length; if trigger wheel is not good
        lbl_adjprf = tk.Label(ProfC_cpane.frame, text= "Adjust Profile", font = HEADING)
        lbl_adjprf.grid(row = 0, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_adjprf = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_adjprf.grid(row = 0, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        AdjProfileButton = tk.Button(ProfC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.adjProfile(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        AdjProfileButton.config(height = 1, width = BTN_GO_W)         
        AdjProfileButton.grid(row=0, column=2, sticky='nsew')
        self.balloon.bind(AdjProfileButton,
                          "Adjust the profile length to \n"
                          "known start and end positions\n"
                          "and/or flip the profile horizontally\n"
                          "(left to right).")

        
        # Set new zero time
        lbl_zt = tk.Label(ProfC_cpane.frame, text= " Set Zero Time", font = HEADING)
        lbl_zt.grid(row = 1, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_zt = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_zt.grid(row = 1, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        SetZeroTimeButton = tk.Button(ProfC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.setZeroTime(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        SetZeroTimeButton.config(height = 1, width = BTN_GO_W)         
        SetZeroTimeButton.grid(row=1, column=2, sticky='nsew')    
        self.balloon.bind(SetZeroTimeButton,
                          "Set the travel time that \n" 
                          "corresponds to the surface.")



        # TimeZero Adjust = align traces
        lbl_altr = tk.Label(ProfC_cpane.frame, text= "Allign Traces", font = HEADING)
        lbl_altr.grid(row = 2, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_altr = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_altr.grid(row = 2, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        TrAlignButton = tk.Button(ProfC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [proj.alignTraces(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        TrAlignButton.config(height = 1, width = BTN_GO_W)         
        TrAlignButton.grid(row=2, column=2, sticky='nsew')
        self.balloon.bind(TrAlignButton,
                         'Automatically shifts each trace up or down\n'
                         'such that the maximum aplitudes of the individual\n'
                         'traces align. Can lead to problems when the maxima\n' 
                         'are not in the air waves. If the results are bad,\n' 
                         'use the "undo" button.')

        

        # truncate Y
        lbl_adjprf = tk.Label(ProfC_cpane.frame, text= "Truncate Y", font = HEADING)
        lbl_adjprf.grid(row = 3, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_adjprf = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_adjprf.grid(row = 3, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        truncYButton = tk.Button(ProfC_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.truncateY(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        truncYButton.config(height = 1, width = BTN_GO_W)         
        truncYButton.grid(row=3, column=2, sticky='nsew')
        self.balloon.bind(truncYButton,
                          "Remove data points at arrival times\n"
                          "later than the chosen value. If velocity\n"
                          "is given: remove data points at depths greater\n"
                          "than the chosen value.")   
        
        #Cut Button
        lbl_cut = tk.Label(ProfC_cpane.frame, text= "Cut Profile", font = HEADING)
        lbl_cut.grid(row = 4, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        
        lbl_cutx = tk.Label(ProfC_cpane.frame, text= "x - value")
        lbl_cutx.grid(row = 5, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_cutx = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_cutx.grid(row = 5, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_cuty = tk.Label(ProfC_cpane.frame, text= "y - value")
        lbl_cuty.grid(row = 5, column = 2, sticky= 'ew', padx= PAD, pady = PAD)
        tb_cuty = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_cuty.grid(row = 5, column = 3, sticky= 'ew', padx= PAD, pady = PAD)

        cutButton = tk.Button(ProfC_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.cut(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        cutButton.config(height = 1, width = BTN_GO_W)         
        cutButton.grid(row=5, column = 4, sticky='nsew')
        self.balloon.bind(cutButton,
                          "trims data to desired along-profile range.") 



        # Dewow
        lbl_dewow = tk.Label(ProfC_cpane.frame, text= "Dewow", font =  HEADING)
        lbl_dewow.grid(row = 6, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_dewow = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_dewow.grid(row = 6, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        DewowButton = tk.Button(ProfC_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.dewow(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        DewowButton.config(height = 1, width = BTN_GO_W)         
        DewowButton.grid(row=6, column=2, sticky='nsew')
        self.balloon.bind(DewowButton,
                          "Trace-wise low-cut filter. Removes\n" 
                          "from each trace a running mean of\n"
                          "chosen window width.")


        # Rem mean trace
        lbl_rmt = tk.Label(ProfC_cpane.frame, text= "Remove Mean Traces", font =  HEADING)
        lbl_rmt.grid(row = 7, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_rmt = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_rmt.grid(row = 7, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        remMeanTraceButton = tk.Button(ProfC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.remMeanTrace(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        remMeanTraceButton.config(height = 1, width = BTN_GO_W)         
        remMeanTraceButton.grid(row=7, column=2, sticky='nsew')
        self.balloon.bind(remMeanTraceButton,
                          "Removes from each trace the average\n" 
                          "of its surrounding traces. This can be\n"
                          "useful to remove air waves, ground\n" 
                          "waves, or horizontal features.")
        
        # Smooth 
        lbl_smth = tk.Label(ProfC_cpane.frame, text= "Smooth (Temp)", font =  HEADING)
        lbl_smth.grid(row = 8, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_smth = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_smth.grid(row = 8, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        SmoothButton = tk.Button(ProfC_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.smooth(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        SmoothButton.config(height = 1, width = BTN_GO_W)         
        SmoothButton.grid(row=8, column=2, sticky='nsew')
        self.balloon.bind(SmoothButton,
                          "Trace-wise high-cut filter. Replaces\n" 
                          "each sample within a trace by a\n"
                          "running mean of chosen window width.")


        # profile Smoothing Button
        lbl_pfsmth = tk.Label(ProfC_cpane.frame, text= "Profile Smoothing", font = HEADING)
        lbl_pfsmth.grid(row = 9, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        
        lbl_trnum = tk.Label(ProfC_cpane.frame, text= "Trace Number")
        lbl_trnum.grid(row = 10, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_trnum = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_trnum.grid(row = 10, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_cpynum = tk.Label(ProfC_cpane.frame, text= "Copies Number")
        lbl_cpynum.grid(row = 10, column = 2, sticky= 'ew', padx= PAD, pady = PAD)
        tb_cpynum = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_cpynum.grid(row = 10, column = 3, sticky= 'ew', padx= PAD, pady = PAD)

        profSmButton = tk.Button(ProfC_cpane.frame,
            text="Go", fg="black",
            command=lambda : [self.profileSmooth(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        profSmButton.config(height = 1, width = BTN_GO_W)         
        profSmButton.grid(row=10, column=4, sticky='nsew')
        self.balloon.bind(profSmButton,
                          "First oversamples the profile (makes 'n' copies\n"
                          "of each trace) and then replaces each trace by\n"
                          "the mean of its neighboring 'm' traces."
                          )
        


        # Gain
        lbl_tpow = tk.Label(ProfC_cpane.frame, text= "Tpow", font =  HEADING)
        lbl_tpow.grid(row = 11, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_smth = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_smth.grid(row = 11, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        tpowButton = tk.Button(ProfC_cpane.frame, 
            text="Go", fg="black",
            command=lambda : [self.tpowGain(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        tpowButton.config(height=1, width=BTN_GO_W)
        tpowButton.grid(row=11, column=2, sticky='nsew')
        self.balloon.bind(tpowButton,
                          "t-power gain. Increases the power of the\n"
                          "signal by a factor of (travel time)^p, where\n"
                          "the user provides p. This gain tends to be\n" 
                          "less aggressive than agc.")

        lbl_agc = tk.Label(ProfC_cpane.frame, text= "AGC", font =  HEADING)
        lbl_agc.grid(row = 12, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_agc = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_agc.grid(row = 12, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        agcButton = tk.Button(
            text="Go",fg="black",
            command=lambda : [self.agcGain(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        agcButton.config(height=1, width=BTN_GO_W)
        agcButton.grid(row=12, column=2, sticky='nsew')
        self.balloon.bind(agcButton,
                          "Automatic gain controll. Normalizes the power\n"
                          "of the signal per given sample window along\n" 
                          "each trace.")

        # show hyperbola
        lbl_hypb = tk.Label(ProfC_cpane.frame, text= "Plot Hyperbola", font = HEADING)
        lbl_hypb.grid(row = 13, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        
        lbl_hypbcent = tk.Label(ProfC_cpane.frame, text= "Hyperbola Center on Profile [m]")
        lbl_hypbcent.grid(row = 14, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_hypbcent = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_hypbcent.grid(row = 14, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_hypbapx = tk.Label(ProfC_cpane.frame, text= "Hyperbola Apex Location [ns]")
        lbl_hypbapx.grid(row = 14, column = 2, sticky= 'ew', padx= PAD, pady = PAD)
        tb_hypbapx = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_hypbapx.grid(row = 14, column = 3, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_hypbvelo = tk.Label(ProfC_cpane.frame, text= "Estimated Velocity [m/ns]")
        lbl_hypbvelo.grid(row = 15, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        tb_hypbvelo = tk.Text(ProfC_cpane.frame, height = 1, width= WIDTH)
        tb_hypbvelo.grid(row = 15, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        hypButton = tk.Button(
            text="show hyperb", fg="black",
            command=lambda : [self.showHyp(proj,a), canvas.draw()])
        hypButton.config(height = 1, width = BTN_GO_W)
        hypButton.grid(row=15, column=2, sticky='nsew')
        self.balloon.bind(hypButton,
                          "Draws a hyperbola depending on profile position,\n"
                          "travel time, and estimated velocity. This can be\n" 
                          "used to find the subsurface velocity when a\n"
                          "a hyperbola is visible in the data. The plotted\n"
                          "hyperbola will disappear when the image is\n" 
                          "refreshed.")



    ####################################################################################### 
    #File operations functions
        
    def loadData(self,proj):
        filename = fd.askopenfilename( filetypes= (("All", "*.*"),
                                                   ("GPRPy (.gpr)", "*.gpr"),
                                                   ("Sensors and Software (.DT1)", "*.DT1"),
                                                   ("GSSI (.DZT)", "*.DZT"),
                                                   ("BSQ header","*.GPRhdr"),
                                                   ("MALA header","*.rad")))
        if filename:
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

    def plotProfileData(self,proj,fig,a,canvas):
        # Clear cursor coordinate cid if if exists to avoid multiple instances
        if 'self.cursor_cid' in locals():
            canvas.mpl_disconnect(self.cursor_cid)            
        dx=proj.profilePos[3]-proj.profilePos[2]
        dt=proj.twtt[3]-proj.twtt[2]
        a.clear()        
        stdcont = np.nanmax(np.abs(proj.data)[:])        
        if proj.velocity is None:
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos)-dx/2.0,
                                                             max(proj.profilePos)+dx/2.0,
                                                             max(proj.twtt)+dt/2.0,
                                                             min(proj.twtt)-dt/2.0],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)
            a.set_ylabel("time [ns]", fontsize=mpl.rcParams['font.size'])
            a.invert_yaxis()
        elif proj.maxTopo is None:
            dy=dt*proj.velocity
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos)-dx/2.0,
                                                             max(proj.profilePos)+dx/2.0,
                                                             max(proj.depth)+dy/2.0,
                                                             min(proj.depth)-dy/2.0],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylabel("depth [m]", fontsize=mpl.rcParams['font.size'])
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)
            a.invert_yaxis()
        else:
            dy=dt*proj.velocity
            a.imshow(proj.data,cmap=self.color.get(),extent=[min(proj.profilePos)-dx/2.0,
                                                             max(proj.profilePos)+dx/2.0,
                                                             proj.minTopo-max(proj.depth)-dy/2.0,
                                                             proj.maxTopo-min(proj.depth)+dy/2.0],
                     aspect="auto",
                     vmin=-stdcont/self.contrast.get(), vmax=stdcont/self.contrast.get())
            a.set_ylabel("elevation [m]", fontsize=mpl.rcParams['font.size'])
            a.set_ylim(self.yrng)
            a.set_xlim(self.xrng)

        a.get_xaxis().set_visible(True)
        a.get_yaxis().set_visible(True)                    
        a.set_xlabel("profile position [m]", fontsize=mpl.rcParams['font.size'])
        a.xaxis.tick_top()
        a.xaxis.set_label_position('top')
        if self.asp is not None:
            a.set_aspect(self.asp)

        # Set grid
        a.grid(self.grid)
            
        # In case you are picking
        figcolsp = 1  # Define the variable figcolsp
        figrowsp = 19 + 1  # Define the variable figrowsp
        if self.picking:
            a.plot(self.picked[:,0],self.picked[:,1],'-x',color='yellow',linewidth=3*self.highfac) 
            a.plot(self.picked[:,0],self.picked[:,1],'-x',color='black',linewidth=2*self.highfac)                               

        # Allow for cursor coordinates being displayed        
        def moved(event):
            if event.xdata is not None and event.ydata is not None:
                canvas.get_tk_widget().itemconfigure(tag, text="(x = %5.5g, y = %5.5g)" % (event.xdata, event.ydata))

        self.cursor_cid = canvas.mpl_connect('button_press_event', moved)
        tag = canvas.get_tk_widget().create_text(20, 20, text="", anchor="nw")

        canvas.get_tk_widget().grid(row=0,column=1,columnspan=figcolsp, rowspan=figrowsp, sticky='nsew')
        canvas.draw()

    def undo(self,proj):
        if self.picking:
            self.picked=self.picked[0:-1,:]
        else:
            proj.undo() 

    def saveData(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".gpr")
        if filename is not '':
            proj.save(filename)

    def exportVTK(self,proj):                    
        outfile = fd.asksaveasfilename()
        if outfile is not '':
            #thickness = sd.askfloat("Input","Profile thickness [m]")
            thickness = 0
            if self.asp is None:
                aspect = 1.0
            else:
                aspect = self.asp
            
            if proj.threeD is None:
                gpyes = mesbox.askyesno("Question","Do you have topography data for this profile?")
                if gpyes:
                    filename = fd.askopenfilename()
                    self.getDelimiter()
                    proj.exportVTK(outfile,gpsinfo=filename,thickness=thickness,delimiter=self.delimiter,aspect=aspect)
            else:
                proj.exportVTK(outfile,gpsinfo=proj.threeD,thickness=thickness,delimiter=self.delimiter,aspect=aspect)
            print('... done with exporting to VTK.')

    #Print Figure
    def printProfileFig(self,proj,fig):
        figname = fd.asksaveasfilename(defaultextension=".pdf")
        if figname is not '':
            dpi = sd.askinteger("Input","Resolution in dots per inch? (Recommended: 600)")
            if dpi is not None:
                fig.savefig(figname, format='pdf', dpi=dpi)        
                # Put what you did in history
                if self.asp is None:
                    histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], dpi=%d)" %(figname,self.color.get(),self.contrast.get(),self.yrng[0],self.yrng[1],self.xrng[0],self.xrng[1],dpi)
                else:
                    histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], asp=%g, dpi=%d)" %(figname,self.color.get(),self.contrast.get(),self.yrng[0],self.yrng[1],self.xrng[0],self.xrng[1],self.asp,dpi)
                proj.history.append(histstr)
        print("Saved figure as %s" %(figname+'.pdf'))

    #For Write Script
    def writeHistory(self,proj):        
        filename = fd.asksaveasfilename(defaultextension=".py")
        if filename is not '':
            proj.writeHistory(filename)
            print("Wrote script to " + filename)

    #######################################################################################
    #View Controls functions
            
    #Full View Function
    def setFullView(self,proj):    
        self.xrng=[np.min(proj.profilePos),np.max(proj.profilePos)]
        if proj.velocity is None:
            self.yrng=[np.min(proj.twtt),np.max(proj.twtt)]
        elif proj.maxTopo is None:
            self.yrng=[np.min(proj.depth),np.max(proj.depth)]
        else:
            self.yrng=[proj.minTopo-np.max(proj.depth),proj.maxTopo-np.min(proj.depth)]
    
    #Grid Function
    def toggleGrid(self):
        self.grid = not self.grid

    #Start Picking Function
    def startPicking(self,proj,fig,a,canvas):
        self.picking = True
        self.picked = np.asmatrix(np.empty((0,2)))
        print("Picking mode on")
        def addPoint(event):
            self.picked = np.append(self.picked,np.asmatrix([event.xdata,event.ydata]),axis=0)
            self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)
            print(self.picked)
        self.pick_cid = canvas.mpl_connect('button_press_event', addPoint)

    #Stop Picking Function
    def stopPicking(self,proj,canvas):
        filename = fd.asksaveasfilename()
        if filename is not '':
            self.picking = False
            canvas.mpl_disconnect(self.pick_cid)
            print("Picking mode off")
            np.savetxt(filename+'_profile.txt',self.picked,delimiter='\t')
            print('saved picked file as "%s"' %(filename+'_profile.txt'))
            # If we have 3D info, also plot it as 3D points
            if proj.threeD is not None:
                # First calculate along-track points
                topoVal = proj.threeD[:,2]
                npos = proj.threeD.shape[0]
                steplen = np.sqrt(
                    np.power( proj.threeD[1:npos,0]-proj.threeD[0:npos-1,0] ,2.0) + 
                    np.power( proj.threeD[1:npos,1]-proj.threeD[0:npos-1,1] ,2.0) +
                    np.power( proj.threeD[1:npos,2]-proj.threeD[0:npos-1,2] ,2.0)
                )
                alongdist = np.cumsum(steplen)
                topoPos = np.append(0,alongdist)
                pick3D = np.zeros((self.picked.shape[0],3))
                # If profile is adjusted, need to start the picked at zero.
                pickProfileShifted = self.picked[:,0] - np.min(proj.profilePos)
                #for i in range(0,3):
                for i in range(0,2):
                    pick3D[:,i] = interp.pchip_interpolate(topoPos,
                                                           proj.threeD[:,i],
                                                           pickProfileShifted).squeeze()
                                                           #self.picked[:,0]).squeeze()
            
                pick3D[:,2] = self.picked[:,1].squeeze()
                    
                np.savetxt(filename+'_3D.txt',pick3D,delimiter='\t')
                print('saved picked file as "%s"' %(filename+'_3D.txt'))       

    #Set X-Range Function
    def setXrng(self):
        xlow = sd.askfloat("Input","Min X value",initialvalue=self.xrng[0])
        if xlow is not None:
            xhigh = sd.askfloat("Input","Max X value",initialvalue=self.xrng[1])
            if xhigh is not None:
                self.xrng=[xlow,xhigh]

    #Set Y-Range Function
    def setYrng(self):
        ylow = sd.askfloat("Input","Min Y value",initialvalue=self.yrng[0])
        if ylow is not None:            
            yhigh = sd.askfloat("Input","Max Y value",initialvalue=self.yrng[1])
            if yhigh is not None:
                self.prevyrng=self.yrng
                self.yrng=[ylow,yhigh]

    #Set Aspect Ratio Function
    def setAspect(self):
        self.asp = sd.askfloat("Input","Plotting aspect ratio", initialvalue=self.asp)

    #######################################################################################
    #Velocity Controls functions
    def setVelocity(self,proj):
        velocity =  sd.askfloat("Input","Radar wave velocity [m/ns]?")        
        if velocity is not None:
            proj.setVelocity(velocity)
            self.prevyrng=self.yrng
            self.yrng=[0,np.max(proj.depth)]

    #Antenna Separation Function
    def antennaSep(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Antenna Sep Error","You have to set the velocity first")
        proj.antennaSep()

    #FK Migration Function
    def fkMigration(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Migration Error","You have to set the velocity first")
        proj.fkMigration()

    #Topo Correct Function
    def topoCorrect(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Topo Correct Error","You have to set the velocity first")
            return
        topofile = fd.askopenfilename()
        if topofile is not '':
            out = self.getDelimiter()    
            proj.topoCorrect(topofile,self.delimiter)
            self.prevyrng=self.yrng
            self.yrng=[proj.minTopo-np.max(proj.depth),proj.maxTopo]
            
    #Profile Smooth Function
    def profileSmooth(self,proj):
        ntraces = sd.askinteger("Input","Smooth over how many traces?")
        if ntraces is not None:
            noversample = sd.askinteger("Input","Make how many copies of each trace?\nRecommended: Same as number of traces to be smoothed.")
            if noversample is not None:
                proj.profileSmooth(ntraces,noversample)

    #Start Picking Function
    def startPicking(self,proj,fig,a,canvas):
        self.picking = True
        self.picked = np.asmatrix(np.empty((0,2)))
        print("Picking mode on")
        def addPoint(event):
            self.picked = np.append(self.picked,np.asmatrix([event.xdata,event.ydata]),axis=0)
            self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)
            print(self.picked)
        self.pick_cid = canvas.mpl_connect('button_press_event', addPoint)

            
    #Stop Picking Function
    def stopPicking(self,proj,canvas):
        filename = fd.asksaveasfilename()
        if filename is not '':
            self.picking = False
            canvas.mpl_disconnect(self.pick_cid)
            print("Picking mode off")
            np.savetxt(filename+'_profile.txt',self.picked,delimiter='\t')
            print('saved picked file as "%s"' %(filename+'_profile.txt'))
            # If we have 3D info, also plot it as 3D points
            if proj.threeD is not None:
                # First calculate along-track points
                topoVal = proj.threeD[:,2]
                npos = proj.threeD.shape[0]
                steplen = np.sqrt(
                    np.power( proj.threeD[1:npos,0]-proj.threeD[0:npos-1,0] ,2.0) + 
                    np.power( proj.threeD[1:npos,1]-proj.threeD[0:npos-1,1] ,2.0) +
                    np.power( proj.threeD[1:npos,2]-proj.threeD[0:npos-1,2] ,2.0)
                )
                alongdist = np.cumsum(steplen)
                topoPos = np.append(0,alongdist)
                pick3D = np.zeros((self.picked.shape[0],3))
                # If profile is adjusted, need to start the picked at zero.
                pickProfileShifted = self.picked[:,0] - np.min(proj.profilePos)
                #for i in range(0,3):
                for i in range(0,2):
                    pick3D[:,i] = interp.pchip_interpolate(topoPos,
                                                           proj.threeD[:,i],
                                                           pickProfileShifted).squeeze()
                                                           #self.picked[:,0]).squeeze()
            
                pick3D[:,2] = self.picked[:,1].squeeze()
                    
                np.savetxt(filename+'_3D.txt',pick3D,delimiter='\t')
                print('saved picked file as "%s"' %(filename+'_3D.txt')) 


    # Show hyperbola
    def showHyp(self,proj,a):
        x0 = sd.askfloat("Input","Hyperbola center on profile [m]", initialvalue=self.hypx)
        if x0 is not None:
            t0 = sd.askfloat("Input","Hyperbola apex location (time [ns])", initialvalue=self.hypt)
            if t0 is not None:
                v  = sd.askfloat("Input","Estimated velocity [m/ns]", initialvalue=self.hypv)
                if v is not None:
                    y=proj.profilePos-x0
                    d=v*t0/2.0
                    k=np.sqrt(d**2 + np.power(y,2))
                    t2=2*k/v
                    a.plot(proj.profilePos,t2,'--c',linewidth=3)
                    self.hypx = x0
                    self.hypt = t0
                    self.hypv = v
    
    # Cut profile
    def cut(self,proj):
        minX = sd.askfloat("Input","Minimum profile position")
        if minX is not None:
            maxX = sd.askfloat("Input","Maximum profile position")
            if maxX is not None:
                proj.cut(minX,maxX)

    def remMeanTrace(self,proj):
        ntraces = sd.askinteger("Input","Remove mean over how many traces?")
        if ntraces is not None:
            proj.remMeanTrace(ntraces=ntraces)


    def tpowGain(self,proj):
        power = sd.askfloat("Input","Power for tpow gain?")
        if power is not None:
            proj.tpowGain(power=power)
        

    def agcGain(self,proj):
        window = sd.askinteger("Input","Window length for AGC?")
        if window is not None:
            proj.agcGain(window=window)

    def truncateY(self,proj):
        maxY = sd.askfloat("Input","Truncate at what y value\n" 
                           "(travel time or depth)")
        if maxY is not None:
            proj.truncateY(maxY)

    def setZeroTime(self,proj):
        newZeroTime = sd.askfloat("Input","New zero time")
        if newZeroTime is not None:
            proj.setZeroTime(newZeroTime=newZeroTime)

    def adjProfile(self,proj):
        flipit = mesbox.askyesno("Question","Flip the profile (left to right)?")
        if flipit:
            proj.flipProfile()        
        minPos = sd.askfloat("Input","Start x coordinate",initialvalue=self.xrng[0])
        if minPos is not None:
            maxPos = sd.askfloat("Input","End x coordinate",initialvalue=self.xrng[1])
            if maxPos is not None:
                proj.adjProfile(minPos=minPos,maxPos=maxPos)
                self.xrng=[minPos,maxPos]
        
    def dewow(self,proj):
        window = sd.askinteger("Input","Dewow window width (number of samples)")
        if window is not None:
            proj.dewow(window=window)


    def smooth(self,proj):
        window = sd.askinteger("Input","Smoothing window width (number of samples)")
        if window is not None:
            proj.smooth(window=window)






    def getDelimiter(self):                
        commaQuery = tk.Toplevel(self.window)
        commaQuery.title("Comma or tab separated?")     
        text = tk.Label(commaQuery,text="Is this a comma- or tab-separated file?",fg='red')
        text.pack(padx=10,pady=10)
        commaButton = tk.Button(commaQuery,text="comma",width=10,
                                command = lambda: [self.setComma(),
                                                   commaQuery.destroy()])
        commaButton.pack(side="left")
        tabButton = tk.Button(commaQuery,text="tab",width=10,
                              command = lambda: [self.setTab(),
                                                 commaQuery.destroy()])
        tabButton.pack(side="right")
        #self.window.frame().tansient(self.window)
        #self.window.frame().grab_set()
        self.window.wait_window(commaQuery)        
    def setComma(self):
        self.delimiter = ','
        print("Delimiter set to comma")
    def setTab(self):
        self.delimiter = '\t'
        print("Delimiter set to tab") 
    


def main():#
	rightcol=9
	figrowsp=19+1
	
	root = tk.Tk()
	GPRPyApp(root)
	root.mainloop()
	
			  
	

if __name__ == '__main__':
	main()
			