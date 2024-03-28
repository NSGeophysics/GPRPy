import sys
import customtkinter as ctk
from tkinter import filedialog, messagebox, ttk
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import simpledialog as sd
from tkinter import messagebox as mesbox
import matplotlib as mpl
from matplotlib import mlab
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
import topLevel_test

import csv
# Importing Collapsible Pane class that we have
# created in separate file
from collapsiblepane import CollapsiblePane as cp
from mayavi import mlab
# import vtkpip
# import matlab.engine

from tvtk.api import tvtk
from vtk.util.numpy_support import vtk_to_numpy
from tkinter import filedialog
import tkinter as tk
 
# Making root window or parent window


# Creating Object of Collapsible Pane Container
# If we do not pass these strings in
# parameter the default strings will appear
# on button that were, expand >>, collapse <<
PAD = 4
WIDTH = 10
HEIGHT = 1
HEADING = ('TkDefaultFont', 13,'bold')
BTN_GO_W = 2


class GPRPyApp:
    

    
    def __init__(self, master):
        MODE = "dark"
        ctk.set_appearance_mode(MODE)
        ctk.set_default_color_theme("dark-blue")
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
        # master.propagate = True


        proj = gp.gprpyProfile()

        btn_frm = ctk.CTkScrollableFrame(master,border_width= 1, width = 570, bg_color="light grey",fg_color="light grey")
    
        btn_frm.grid(row = 0, column= 0, sticky= 'ns', rowspan = 2)
        # btn_frm.grid_propagate(0)
        btn_frm.columnconfigure(1, minsize=570, weight=1)
        btn_frm.grid_propagate = 1
        master.rowconfigure(1, weight=10) 

        fig=Figure(figsize=(self.widfac,self.highfac),facecolor="#d3d3d3")
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
        FM_cpane.grid_propagate = True

        #Set frame color:
        # Initialize style
        s = ttk.Style()
        # Create style used by default for all Frames
        s.configure('TFrame', background='#d3d3d3')
        s.configure('TButton', font=('Helvetica', 12))
        


        LoadButton = ctk.CTkButton(FM_cpane.frame,
            text="Import Data", 
            command=lambda : [self.loadData(proj), 
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        LoadButton.configure(height = HEIGHT, width = WIDTH)         
        LoadButton.grid(row=0, column=0, sticky='nsew',pady = PAD)
        self.balloon.bind(LoadButton,"Load .gpr, .DT1, or .DZT data.")


        undoButton = ctk.CTkButton(FM_cpane.frame,
            text="Undo",
            command=lambda : [self.resetYrng(proj),
                              self.undo(proj), 
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        undoButton.configure(height = HEIGHT, width = WIDTH)
        undoButton.grid(row=1, column=0, sticky='ew',pady= PAD )
        self.balloon.bind(undoButton,
                          '"Undoes" the most recent processing step and\n' 
                          'sets the data back to its previous state.\n' 
                          'This also removes the most recent processing\n'
                          'step from the history. Does not revert\n' 
                          'visualization settings such as "set x-range"\n'
                          'etc.')

        

        SaveButton = ctk.CTkButton(FM_cpane.frame,
            text="Save Data", 
            command=lambda : self.saveData(proj))
        SaveButton.configure(height = HEIGHT, width = WIDTH)         
        SaveButton.grid(row= 4, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(SaveButton,
                          'saves the processed data including its history in a\n'
                          '.gpr file. The resulting file will contain absolute\n'
                          'path names of the used data and topography files.\n'
                          'Visualization settings such as "set x-range" or\n'
                          '"contrast" will not be saved.')

        PrintButton = ctk.CTkButton(FM_cpane.frame,
            text="Print Figure", 
            command=lambda : self.printProfileFig(proj=proj,fig=fig))
        PrintButton.configure(height = HEIGHT, width = WIDTH)         
        PrintButton.grid(row=5, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(PrintButton,
                          "Saves the current visible figure in a pdf with \n"
                          "chosen resolution. If there is a hyperbola on\n" 
                          "the current figure, then the hyperbola will also\n"
                          "appear on the printed figure.")


        # Export to VTK
        VTKButton = ctk.CTkButton(FM_cpane.frame,
            text="Export to VTK",
            command = lambda : self.exportVTK(proj))
        VTKButton.configure(height = HEIGHT, width = WIDTH)
        VTKButton.grid(row=6, column=0, sticky='ew',pady = PAD)
        self.balloon.bind(VTKButton,
                          "Exports the processed figure to a\n"
                          "VTK format, that can be read by\n" 
                          "Paraview or similar 3D programs.")
        


        
        # Write script
        HistButton = ctk.CTkButton(FM_cpane.frame,
            text="Write Script",
            command=lambda : self.writeHistory(proj))
        HistButton.configure(height = HEIGHT, width = WIDTH)         
        HistButton.grid(row=7, column=0, sticky='ew',pady = PAD)
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
        VC_cpane.grid_propagate = True
        VC_cpane.columnconfigure(1, minsize=570, weight=1)



        # Full view
        FullButton = ctk.CTkButton(VC_cpane.frame,
            text="Full View", 
            command=lambda : [self.setFullView(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        FullButton.configure(height = HEIGHT, width = WIDTH)         
        FullButton.grid(row=0, column=0, sticky='ew', padx = PAD, pady = PAD)
        self.balloon.bind(FullButton,"Resets x- and y-axis limits to full data.")

        
        # Grid button
        GridButton = ctk.CTkButton(VC_cpane.frame,
            text="Grid", 
            command=lambda : [self.toggleGrid(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        GridButton.configure(height = HEIGHT, width = WIDTH)         
        GridButton.grid(row = 0, column=1, sticky='ew', padx = PAD, pady = PAD)
        self.balloon.bind(GridButton,"Toggles grid on/off.")

        startPickButton = ctk.CTkButton(VC_cpane.frame,
            text="start pick", 
            command=lambda : self.startPicking(proj,fig=fig,a=a,canvas=canvas))        
        startPickButton.configure(height = HEIGHT, width = WIDTH) 
        startPickButton.grid(row=1, column=0, sticky='ew',padx = PAD, pady = PAD)
        self.balloon.bind(startPickButton,
                          "Start collecting location information\n" 
                          "by clicking on the profile.")  
        

        stopPickButton = ctk.CTkButton(VC_cpane.frame,
            text="stop pick", 
            command=lambda : [self.stopPicking(proj,canvas),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        stopPickButton.configure(height = HEIGHT, width = WIDTH) 
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

        #Changed to take input from text box
        self.tb_minx = tk.Text(VC_cpane.frame, height = 1, width = 3)
        self.tb_minx.grid(row = 3, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_maxx = tk.Label(VC_cpane.frame, text= "Max:")
        lbl_maxx.grid(row = 3, column = 2, sticky= 'ew', padx= PAD, pady = PAD)

        #Changed to take input from text box
        self.tb_maxx = tk.Text(VC_cpane.frame, height = 1, width = 3)
        self.tb_maxx.grid(row = 3, column = 3, sticky= 'ew', padx= PAD, pady = PAD)


        XrngButton = ctk.CTkButton(VC_cpane.frame,
            text="Go", 
            command=lambda : [self.setXrng(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        XrngButton.configure(height = HEIGHT, width = 3)      
        XrngButton.grid(row=3, column=4, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(XrngButton,"Set the x-axis display limits.")
        

        # Y range
        lbl_yrng = tk.Label(VC_cpane.frame, text= "Set Y - Range", font = HEADING)
        lbl_yrng.grid(row = 4, column = 0, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_miny = tk.Label(VC_cpane.frame, text= "Min:")
        lbl_miny.grid(row = 5, column = 0, sticky= 'ew', padx= PAD, pady = PAD)

        self.tb_miny = tk.Text(VC_cpane.frame, height = 1, width = 3)
        self.tb_miny.grid(row = 5, column = 1, sticky= 'ew', padx= PAD, pady = PAD)

        lbl_maxy = tk.Label(VC_cpane.frame, text= "Max:")
        lbl_maxy.grid(row = 5, column = 2, sticky= 'ew', padx= PAD, pady = PAD)

        self.tb_maxy = tk.Text(VC_cpane.frame, height = 1, width = 3)
        self.tb_maxy.grid(row = 5, column = 3, sticky= 'ew', padx= PAD, pady = PAD)


        YrngButton = ctk.CTkButton(VC_cpane.frame,
            text="Go", 
            command=lambda : [self.setYrng(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        YrngButton.configure(height = HEIGHT, width = 3)         
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
        
        # Go button for colour and contrast
        PlotButton = ctk.CTkButton(VC_cpane.frame,
            text="Go", 
            command=lambda : [self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        PlotButton.configure(height = HEIGHT, width = 3)         
        PlotButton.grid(row=6, column= 3, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(PlotButton,"Set Colour and Contrast")
    
        # Aspect
        lbl_aspr = tk.Label(VC_cpane.frame, text= "Aspect Ratio", font = HEADING)
        lbl_aspr.grid(row = 7, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        self.tb_aspr = tk.Text(VC_cpane.frame, height=1, width=3)
        self.tb_aspr.grid(row=7, column=1, sticky='ew', padx=PAD, pady=PAD)
        AspButton = ctk.CTkButton(VC_cpane.frame,
            text="Go", 
            command=lambda : [self.setAspect(),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])                              
        AspButton.configure(height = HEIGHT, width = 3)      
        AspButton.grid(row=7, column=3, sticky='nsew',padx = PAD, pady = PAD)
        self.balloon.bind(lbl_aspr, "Set the aspect ratio between x- and y-axis.")
    
###############################################################################################
        #View Control Pane and all buttons within it
        Velo_cpane = cp(btn_frm, 'Velocity Controls -', 'Velocity Controls +')
        Velo_cpane.grid(row = 2, column = 0, sticky = 'ew',) 

        # Set Velocity
        lbl_velo = tk.Label(Velo_cpane.frame, text= "Set Velocity", font = HEADING)
        lbl_velo.grid(row = 0, column = 0, sticky= 'ew', padx= PAD, pady = PAD)
        self.tb_velo = tk.Text(Velo_cpane.frame, height = 1, width= WIDTH)
        self.tb_velo.grid(row = 0, column = 1, sticky= 'ew', padx= PAD, pady = PAD)
        setVelButton = ctk.CTkButton(Velo_cpane.frame, 
            text="Go", 
            command=lambda : [self.setVelocity(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        setVelButton.configure(height = 1, width = 3)         
        setVelButton.grid(row=0, column=3, sticky='nsew')
        self.balloon.bind(setVelButton,
                          "Set the known subsurface radar velocity. This\n" 
                          "turns the y-axis from travel time to depth.\n"
                          "This step is necessary for topographic correction.")


        
        # Correct for antenna separation
        antennaSepButton = ctk.CTkButton(Velo_cpane.frame,
            text="Antenna Sep Correct", 
            command=lambda : [self.antennaSep(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        antennaSepButton.configure(height = HEIGHT, width = WIDTH+10)         
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
        migButton = ctk.CTkButton(Velo_cpane.frame, 
            text="Go", 
            command=lambda : [self.fkMigration(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        migButton.configure(height = 1, width = 3)         
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
        topoCorrectButton = ctk.CTkButton(Velo_cpane.frame,
            text="Go", 
            command=lambda : [self.topoCorrect(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        topoCorrectButton.configure(height = 1, width = 3)
        topoCorrectButton.grid(row=3, column=3, sticky='nsew')
        self.balloon.bind(topoCorrectButton,
                          "Reads a comma- or tab-separated file containing\n" 
                          "either 3 columns (easting, northing, elevation)\n" 
                          "or two columns (profile position, elevation).\n" 
                          "All coordinates in meters.")       


####################################################################################### 
        # Profile Controls
        ProfC_cpane = cp(btn_frm, 'Profile Controls -', 'Profile Controls +')
        ProfC_cpane.grid(row=3, column=0, sticky='ew')
        ProfC_cpane.columnconfigure(2,weight=100)



        # Adjust Profile
        lbl_adjprf = tk.Label(ProfC_cpane.frame, text="Adjust Profile", font=HEADING)
        lbl_adjprf.grid(row=0, column=0, sticky='ew', padx=PAD, pady=PAD)

        lbl_adjprfMin = tk.Label(ProfC_cpane.frame, text="Min X-Value:")
        lbl_adjprfMin.grid(row=1, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_adjprfMin = tk.Text(ProfC_cpane.frame, height=1, width=3)  # Reduced width to 2
        self.tb_adjprfMin.grid(row=1, column=1, sticky='ew', padx=PAD, pady=PAD)

        lbl_adjprfMax = tk.Label(ProfC_cpane.frame, text="Max X- Value:")
        lbl_adjprfMax.grid(row=1, column=2, sticky='ew', padx=PAD, pady=PAD)
        self.tb_adjprfMax = tk.Text(ProfC_cpane.frame, height=1, width=3)  # Kept width as 3
        self.tb_adjprfMax.grid(row=1, column=3, sticky='ew', padx=PAD, pady=PAD)

        AdjProfileButton = ctk.CTkButton(ProfC_cpane.frame,
            text="Go", 
            command=lambda : [self.adjProfile(proj),
                              self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)])
        AdjProfileButton.configure(height = 1, width = 2)         
        AdjProfileButton.grid(row=1, column=4, sticky='nsew')
        self.balloon.bind(AdjProfileButton,
                          "Adjust the profile length to \n"
                          "known start and end positions\n"
                          "and/or flip the profile horizontally\n"
                          "(left to right).")

        
        # Set new zero time
        lbl_zt = tk.Label(ProfC_cpane.frame, text=" Set Zero Time", font=HEADING)
        lbl_zt.grid(row=2, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_zt = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_zt.grid(row=2, column=1, sticky='ew', padx=PAD, pady=PAD)

        SetZeroTimeButton = ctk.CTkButton(ProfC_cpane.frame,
                                    text="Go",
                                    command=lambda: [self.setZeroTime(proj),
                                                    self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        SetZeroTimeButton.configure(height=1, width=BTN_GO_W)
        SetZeroTimeButton.grid(row=2, column=2, sticky='nsew')
        self.balloon.bind(SetZeroTimeButton,
                        "Set the travel time that \n"
                        "corresponds to the surface.")


        # TimeZero Adjust = align traces
        lbl_altr = tk.Label(ProfC_cpane.frame, text="Allign Traces", font=HEADING)
        lbl_altr.grid(row=3, column=0, sticky='ew', padx=PAD, pady=PAD)
        tb_altr = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        tb_altr.grid(row=3, column=1, sticky='ew', padx=PAD, pady=PAD)

        TrAlignButton = ctk.CTkButton(ProfC_cpane.frame,
                                text="Go", 
                                command=lambda: [proj.alignTraces(),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        TrAlignButton.configure(height=1, width=BTN_GO_W)
        TrAlignButton.grid(row=3, column=2, sticky='nsew')
        self.balloon.bind(TrAlignButton,
                        'Automatically shifts each trace up or down\n'
                        'such that the maximum aplitudes of the individual\n'
                        'traces align. Can lead to problems when the maxima\n'
                        'are not in the air waves. If the results are bad,\n'
                        'use the "undo" button.')


        # truncate Y
        lbl_adjprf = tk.Label(ProfC_cpane.frame, text="Truncate Y", font=HEADING)
        lbl_adjprf.grid(row=4, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_adjprf = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_adjprf.grid(row=4, column=1, sticky='ew', padx=PAD, pady=PAD)

        truncYButton = ctk.CTkButton(ProfC_cpane.frame,
                                text="Go", 
                                command=lambda: [self.truncateY(proj),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        truncYButton.configure(height=1, width=BTN_GO_W)
        truncYButton.grid(row=4, column=2, sticky='nsew')
        self.balloon.bind(truncYButton,
                        "Remove data points at arrival times\n"
                        "later than the chosen value. If velocity\n"
                        "is given: remove data points at depths greater\n"
                        "than the chosen value.")

        # Cut Button
        lbl_cut = tk.Label(ProfC_cpane.frame, text="Cut Profile", font=HEADING)
        lbl_cut.grid(row=5, column=0, sticky='ew', padx=PAD, pady=PAD)

        lbl_cutx = tk.Label(ProfC_cpane.frame, text="x - value")
        lbl_cutx.grid(row=6, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_cutx = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_cutx.grid(row=6, column=1, sticky='ew', padx=PAD, pady=PAD)

        lbl_cuty = tk.Label(ProfC_cpane.frame, text="y - value")
        lbl_cuty.grid(row=6, column=2, sticky='ew', padx=PAD, pady=PAD)
        self.tb_cuty = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_cuty.grid(row=6, column=3, sticky='ew', padx=PAD, pady=PAD)

        cutButton = ctk.CTkButton(ProfC_cpane.frame,
                            text="Go", 
                            command=lambda: [self.cut(proj),
                                            self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        cutButton.configure(height=1, width=BTN_GO_W)
        cutButton.grid(row=6, column=4, sticky='nsew')
        self.balloon.bind(cutButton,
                        "trims data to desired along-profile range.")


        # Dewow
        lbl_dewow = tk.Label(ProfC_cpane.frame, text="Dewow", font=HEADING)
        lbl_dewow.grid(row=7, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_dewow = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_dewow.grid(row=7, column=1, sticky='ew', padx=PAD, pady=PAD)

        DewowButton = ctk.CTkButton(ProfC_cpane.frame,
                                text="Go", 
                                command=lambda: [self.dewow(proj),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        DewowButton.configure(height=1, width=BTN_GO_W)
        DewowButton.grid(row=7, column=2, sticky='nsew')
        self.balloon.bind(DewowButton,
                        "Trace-wise low-cut filter. Removes\n"
                        "from each trace a running mean of\n"
                        "chosen window width.")


        # Rem mean trace
        lbl_rmt = tk.Label(ProfC_cpane.frame, text="Remove Mean Traces", font=HEADING)
        lbl_rmt.grid(row=8, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_rmt = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_rmt.grid(row=8, column=1, sticky='ew', padx=PAD, pady=PAD)

        remMeanTraceButton = ctk.CTkButton(ProfC_cpane.frame,
                                    text="Go", 
                                    command=lambda: [self.remMeanTrace(proj),
                                                        self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        remMeanTraceButton.configure(height=1, width=BTN_GO_W)
        remMeanTraceButton.grid(row=8, column=2, sticky='nsew')
        self.balloon.bind(remMeanTraceButton,
                        "Removes from each trace the average\n"
                        "of its surrounding traces. This can be\n"
                        "useful to remove air waves, ground\n"
                        "waves, or horizontal features.")

        # Smooth
        lbl_smth = tk.Label(ProfC_cpane.frame, text="Smooth (Temp)", font=HEADING)
        lbl_smth.grid(row=9, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_smth1 = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_smth1.grid(row=9, column=1, sticky='ew', padx=PAD, pady=PAD)

        SmoothButton = ctk.CTkButton(ProfC_cpane.frame,
                                text="Go", 
                                command=lambda: [self.smooth(proj),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        SmoothButton.configure(height=1, width=BTN_GO_W)
        SmoothButton.grid(row=9, column=2, sticky='nsew')
        self.balloon.bind(SmoothButton,
                        "Trace-wise high-cut filter. Replaces\n"
                        "each sample within a trace by a\n"
                        "running mean of chosen window width.")


        # profile Smoothing Button
        lbl_pfsmth = tk.Label(ProfC_cpane.frame, text="Profile Smoothing", font=HEADING)
        lbl_pfsmth.grid(row=10, column=0, sticky='ew', padx=PAD, pady=PAD)

        lbl_trnum = tk.Label(ProfC_cpane.frame, text="Trace Number")
        lbl_trnum.grid(row=11, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_trnum = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_trnum.grid(row=11, column=1, sticky='ew', padx=PAD, pady=PAD)

        lbl_cpynum = tk.Label(ProfC_cpane.frame, text="Copies Number")
        lbl_cpynum.grid(row=11, column=2, sticky='ew', padx=PAD, pady=PAD)
        self.tb_cpynum = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_cpynum.grid(row=11, column=3, sticky='ew', padx=PAD, pady=PAD)

        profSmButton = ctk.CTkButton(ProfC_cpane.frame,
                                text="Go", 
                                command=lambda: [self.profileSmooth(proj),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        profSmButton.configure(height=1, width=BTN_GO_W)
        profSmButton.grid(row=11, column=4, sticky='nsew')
        self.balloon.bind(profSmButton,
                        "First oversamples the profile (makes 'n' copies\n"
                        "of each trace) and then replaces each trace by\n"
                        "the mean of its neighboring 'm' traces."
                        )


        # Gain
        lbl_tpow = tk.Label(ProfC_cpane.frame, text="Tpow", font=HEADING)
        lbl_tpow.grid(row=12, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_smth = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_smth.grid(row=12, column=1, sticky='ew', padx=PAD, pady=PAD)

        tpowButton = ctk.CTkButton(ProfC_cpane.frame,
                            text="Go", 
                            command=lambda: [self.tpowGain(proj),
                                                self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        tpowButton.configure(height=1, width=BTN_GO_W)
        tpowButton.grid(row=12, column=2, sticky='nsew')
        self.balloon.bind(tpowButton,
                        "t-power gain. Increases the power of the\n"
                        "signal by a factor of (travel time)^p, where\n"
                        "the user provides p. This gain tends to be\n"
                        "less aggressive than agc.")

        lbl_agc = tk.Label(ProfC_cpane.frame, text="AGC", font=HEADING)
        lbl_agc.grid(row=13, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_agc = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_agc.grid(row=13, column=1, sticky='ew', padx=PAD, pady=PAD)

        agcButton = ctk.CTkButton(ProfC_cpane.frame,
                            text="Go", 
                            command=lambda: [self.agcGain(proj),
                                            self.plotProfileData(proj, fig=fig, a=a, canvas=canvas)])
        agcButton.configure(height=1, width=BTN_GO_W)
        agcButton.grid(row=13, column=2, sticky='nsew')
        self.balloon.bind(agcButton,
                        "Automatic gain controll. Normalizes the power\n"
                        "of the signal per given sample window along\n"
                        "each trace.")

        # show hyperbola
        lbl_hypb = tk.Label(ProfC_cpane.frame, text="Plot Hyperbola", font=HEADING)
        lbl_hypb.grid(row=14, column=0, sticky='ew', padx=PAD, pady=PAD)

        lbl_hypbcent = tk.Label(ProfC_cpane.frame, text="Hyperbola Center on Profile [m]")
        lbl_hypbcent.grid(row=15, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_hypbcent = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_hypbcent.grid(row=15, column=1, sticky='ew', padx=PAD, pady=PAD)

        lbl_hypbapx = tk.Label(ProfC_cpane.frame, text="Hyperbola Apex Location [ns]")
        lbl_hypbapx.grid(row=15, column=2, sticky='ew', padx=PAD, pady=PAD)
        self.tb_hypbapx = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_hypbapx.grid(row=15, column=3, sticky='ew', padx=PAD, pady=PAD)

        lbl_hypbvelo = tk.Label(ProfC_cpane.frame, text="Estimated Velocity [m/ns]")
        lbl_hypbvelo.grid(row=16, column=0, sticky='ew', padx=PAD, pady=PAD)
        self.tb_hypbvelo = tk.Text(ProfC_cpane.frame, height=1, width=WIDTH)
        self.tb_hypbvelo.grid(row=16, column=1, sticky='ew', padx=PAD, pady=PAD)

        hypButton = ctk.CTkButton(ProfC_cpane.frame,
                            text="show hyperb", 
                            command=lambda: [self.showHyp(proj, a), canvas.draw()])
        hypButton.configure(height=1, width=BTN_GO_W)
        hypButton.grid(row=16, column=2, sticky='nsew')
        self.balloon.bind(hypButton,
                        "Draws a hyperbola depending on profile position,\n"
                        "travel time, and estimated velocity. This can be\n"
                        "used to find the subsurface velocity when a\n"
                        "a hyperbola is visible in the data. The plotted\n"
                        "hyperbola will disappear when the image is\n"
                        "refreshed.")
        
        # # Interpolation CP =====================================================================
        Interp_cpane = cp(btn_frm, 'Interpolate -', 'Interpolate +')
        Interp_cpane.grid(row=5, column=0, sticky='ew')

        makeDataCubeBtn = ctk.CTkButton(Interp_cpane.frame,
            text="Interpolate",
            command=self.loadInfoCollectScreen) 
        makeDataCubeBtn.configure(height=HEIGHT, width=WIDTH)         
        makeDataCubeBtn.grid(row=0, column=0, sticky='nsew', pady=PAD)
        self.balloon.bind(makeDataCubeBtn, "Click to open the interface to create datacube.")



        # DataCube Button CP ==========================================================================
        DataC_cpane = cp(btn_frm, 'DataCube -', 'DataCube +')
        DataC_cpane.grid(row=4, column=0, sticky='ew')
        DataC_cpane.grid_propagate
        DataC_cpane.pack_propagate

        importDataCubeBtn = ctk.CTkButton(DataC_cpane.frame,
            text="Import VTS/VTK File", 
            command=self.load_vtk_file) 
        importDataCubeBtn.configure(height=HEIGHT, width=WIDTH)         
        importDataCubeBtn.grid(row=0, column=0, sticky='nsew', pady=PAD)

        # Inside your GUI setup where you define buttons
        generate_data_btn = ctk.CTkButton(DataC_cpane.frame, text="Generate Data", command=self.generate_data_for_current_position)
        generate_data_btn.configure(height=HEIGHT, width=WIDTH)
        generate_data_btn.grid(row=1, column=0, sticky='nsew', pady=PAD)

        # Reset Data Points list for new slice
        new_slice_button = ctk.CTkButton(DataC_cpane.frame, 
                                     text="Start New Slice",  command=self.reset_data_points_for_new_slice)
        new_slice_button.configure(height=HEIGHT, width=WIDTH)
        new_slice_button.grid(row=2, column=0, sticky='nsew', pady=PAD)
        
        # Export Data Points to CSV
        exportBtn = ctk.CTkButton(DataC_cpane.frame,
                              text="ExportCsv", 
                              command=self.prompt_export_data_points) 
        exportBtn.configure(height=HEIGHT, width=WIDTH)
        exportBtn.grid(row=3, column=0, sticky='nsew', pady=PAD)
        

        
        self.data_points = []
        self.scalar_cut_plane = None 

    def loadInfoCollectScreen(self):
        ld = topLevel_test.Info_Collect()


    def generate_data_for_current_position(self):
        if not hasattr(self, 'scalar_cut_plane'):
            print("Scalar cut plane is not initialized.")
            return

        current_position = self.scalar_cut_plane.implicit_plane.origin
        print(f"Generating data for point:{ current_position}")

        z_depth = current_position[2]  # Use the Z component
        self.generate_slice_data(z_depth)
        print(f"Data generated for Z={z_depth} with {len(self.data_points)} points.")

    def generate_slice_data(self, z_depth):
        x_start, x_end, x_interval = 0, 10, 0.02  # Define grid boundaries and step size
        y_start, y_end, y_interval = 0, 10, 0.02

        self.data_points = []  # Clear existing points for a new slice

        for x in np.arange(x_start, x_end, x_interval):
            for y in np.arange(y_start, y_end, y_interval):
                if len(self.data_points) >= 500:
                    print(f"Reached 500 data points limit. Stopping data generation.")
                    return  # Stop adding more points once 500 have been added
                g_value = self.get_scalar_value_at_point(x, y, z_depth)
                if g_value is not None:
                    self.data_points.append((x, y, z_depth, g_value))
                    # Temporary print statement for debugging:
                    print(f"X: {x}, Y: {y}, Z: {z_depth}, G: {g_value}")

        print(f"Generated slice at Z={z_depth} with {len(self.data_points)} points.")
        
    def reset_data_points_for_new_slice(self):
        self.data_points = []  # Reset the list to empty
        print("Data points reset for a new slice.")
        
    def load_vtk_file(self):
        # Open a file dialog to select the VTK file
        self.vtk_file = filedialog.askopenfilename(filetypes=[("VTK files", "*.vtk"), ("VTS files", "*.vts")])
        if self.vtk_file:
            # Call function to display the VTK file
            self.display_vtk_file()


    def move_cut_plane(self, direction, step=2.0):
        if not hasattr(self, 'scalar_cut_plane'):
            print("Scalar cut plane is not initialized.")
            return

        # Calculate the new Z depth based on direction
        current_position = self.scalar_cut_plane.implicit_plane.origin
        if direction == 'up':
            new_z = current_position[2] + step / 100.0  # Convert cm to meters if necessary
        elif direction == 'down':
            new_z = current_position[2] - step / 100.0
        else:
            print("Invalid direction. Use 'up' or 'down'.")
            return

        # Update the plane's Z position
        self.scalar_cut_plane.implicit_plane.origin = (current_position[0], current_position[1], new_z)

        # Create a scalar cut plane
        # scalar_cut_plane = mlab.pipeline.scalar_cut_plane(data)    
        # Force update visualization
        self.scalar_cut_plane.implicit_plane.scene.render()

        # Generate data points for the new slice
        #self.generate_slice_data(new_z)

    
                
    def display_vtk_file(self):
        src = mlab.pipeline.open(self.vtk_file)

        try:
            if hasattr(src.outputs[0], 'input'):
                vtk_dataset = src.outputs[0].input  # Accessing the input to AssignAttribute
                if hasattr(vtk_dataset, 'point_data'):
                    scalar_range = vtk_dataset.point_data.scalars.range
                else:
                    raise AttributeError("Failed to find point_data in vtk_dataset.")
            else:
                raise AttributeError("src.outputs[0] does not have 'input'.")
        except AttributeError as e:
            print(e)
            return  # Exit if we cannot correctly access the dataset

        # Generate contour levels within the valid scalar range of the data
        contour_levels = [scalar_range[0] + (scalar_range[1] - scalar_range[0]) * 0.25,
                        scalar_range[0] + (scalar_range[1] - scalar_range[0]) * 0.75]

        mlab.pipeline.iso_surface(src, contours=contour_levels, opacity=0.3)
        self.scalar_cut_plane = mlab.pipeline.scalar_cut_plane(src, plane_orientation='z_axes')

        # mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', nb_labels=5)
        # # Retrieve the bounds of your dataset for grid alignment

        # mlab.draw()

        # After initializing self.scalar_cut_plane
        self.scalar_cut_plane.implicit_plane.widget.add_observer('InteractionEvent', self.on_cut_plane_move)

        # Inside display_vtk_file, after successfully loading the dataset
        self.current_vtk_dataset = vtk_dataset  # Ensure this references your loaded VTK dataset object


    def on_cut_plane_move(self, obj, evt):
        origin = self.scalar_cut_plane.implicit_plane.origin
        # Decompose the origin into x, y, and z components
        x, y, z = origin

        # Call get_scalar_value_at_point with x, y, z as separate arguments
        g_value = self.get_scalar_value_at_point(x, y, z)
        
        # Store the coordinates and scalar value
        self.data_points.append((x, y, z, g_value))
        
        # Print the current X, Y, Z, and G values
        print(f"X: {x}, Y: {y}, Z: {z}, G: {g_value}")

        # Optionally update text in the visualization to show current coordinates
        # This part assumes you have a method or logic to manage the text display correctly
        # Make sure to check or implement logic to avoid overlaying multiple text instances
        self.update_coord_text(x, y, z, g_value)

    def update_coord_text(self, x, y, z, g_value):
        # Assuming 'coord_text' is an attribute to store the text object
        # Check if it exists and remove it before creating a new one to avoid overlay
        if hasattr(self, 'coord_text'):
            self.coord_text.remove()

        # Display new coordinates as text
        self.coord_text = mlab.text(0.01, 0.01, f"X: {x:.2f}, Y: {y:.2f}, Z: {z:.2f}, G: {g_value:.2f}", width=0.25)
        mlab.draw()


    def get_scalar_value_at_point(self, x, y, z):
        # Ensure 'self.current_vtk_dataset' is properly set to your VTK dataset
        vtk_data = self.current_vtk_dataset

        # Assuming 'vtk_data' has been properly set up to refer to your current dataset
        point = [x, y, z]
        point_data = tvtk.PolyData(points=[point])
        probe = tvtk.ProbeFilter()
        probe.set_input_data(point_data)
        probe.set_source_data(vtk_data)
        probe.update()

        if probe.output.point_data.scalars:
            scalar_value = probe.output.point_data.scalars.to_array()[0]
            return scalar_value
        else:
            return None  # No scalar value was interpolated at this point

        
    def automated_data_sampling(self):
        x_start, x_end, x_interval = 0, 10, 0.02
        y_start, y_end, y_interval = 0, 10, 0.02
        z_depth = 0.5  # Example depth, adjust based on your requirements

        self.reset_data_points_for_new_slice()

        for x in np.arange(x_start, x_end, x_interval):
            for y in np.arange(y_start, y_end, y_interval):
                g_value = self.get_scalar_value_at_point(x, y, z_depth)
                if g_value is not None:
                    self.data_points.append((x, y, z_depth, g_value))

        print(f"Collected {len(self.data_points)} data points for the new slice.")
        

    def get_scalar_value_at_point(self, x, y, z):
        # This implementation needs to access the scalar value from your dataset at the (x, y, z) location.
        # The actual implementation will depend on how your dataset is structured.
        # The following is a placeholder that assumes you can directly query your dataset.

        vtk_data = self.current_vtk_dataset  # Ensure this is set to your dataset when loaded.

        # Construct a point and use a probe to get the scalar value at this point
        point = [x, y, z]
        point_data = tvtk.PolyData(points=[point])
        probe = tvtk.ProbeFilter()
        probe.set_input_data(point_data)
        probe.set_source_data(vtk_data)
        probe.update()

        if probe.output.point_data.scalars:
            scalar_value = probe.output.point_data.scalars.to_array()[0]
            return scalar_value
        else:
            return None  # Handle cases where no data is interpolated at the point


    def prompt_export_data_points(self):
        # Prompt the user for a file path to export
        filepath = filedialog.asksaveasfilename(defaultextension=".csv",
                                                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if filepath:
            self.export_data_points(filepath)


    def export_data_points(self, filepath):
            with open(filepath, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['X', 'Y', 'Z', 'G'])  # Header
                for x, y, z, g_value in self.data_points:
                    # Format x, y, z, and g_value to have at least 4 significant digits
                    writer.writerow([f"{x:.4g}", f"{y:.4g}", f"{z:.16g}", f"{g_value:.16g}"])
            print(f"Data points exported to {filepath}")

    # if the data is a strucutred grid we use this 
    # def get_scalar_value_at_point(self, src, point):
    #     # Assuming src is a structured grid source for simplicity
    #     vtk_data = src.outputs[0]
    #     if isinstance(vtk_data, tvtk.StructuredGrid):
    #         x, y, z = point
    #         spacing = vtk_data.spacing
    #         origin = vtk_data.origin
    #         i = int(round((x - origin[0]) / spacing[0]))
    #         j = int(round((y - origin[1]) / spacing[1]))
    #         k = int(round((z - origin[2]) / spacing[2]))
    #         try:
    #             scalar_value = vtk_data.point_data.scalars[i, j, k]
    #             return scalar_value
    #         except IndexError:
    #             return None  # Point is outside the grid
    #     return None  # Placeholder for non-structured grid data

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
        if filename != '':
            proj.save(filename)

    def exportVTK(self,proj):                    
        outfile = fd.asksaveasfilename()
        if outfile != '':
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
        if figname != '':
            dpi = sd.askinteger("Input","Resolution in dots per inch? (Recommended: 600)")
            if dpi != None:
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
        if filename != '':
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
        if filename != '':
            self.picking = False
            canvas.mpl_disconnect(self.pick_cid)
            print("Picking mode off")
            np.savetxt(filename+'_profile.txt',self.picked,delimiter='\t')
            print('saved picked file as "%s"' %(filename+'_profile.txt'))
            # If we have 3D info, also plot it as 3D points
            if proj.threeD != None:
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
    # def setXrng(self):
    #     xlow = sd.askfloat("Input","Min X value",initialvalue=self.xrng[0])
    #     if xlow is not None:
    #         xhigh = sd.askfloat("Input","Max X value",initialvalue=self.xrng[1])
    #         if xhigh is not None:
    #             self.xrng=[xlow,xhigh]

    # #Set Y-Range Function
    # def setYrng(self):
    #     ylow = sd.askfloat("Input","Min Y value",initialvalue=self.yrng[0])
    #     if ylow is not None:            
    #         yhigh = sd.askfloat("Input","Max Y value",initialvalue=self.yrng[1])
    #         if yhigh is not None:
    #             self.prevyrng=self.yrng
    #             self.yrng=[ylow,yhigh]

    def setXrng(self):
        xlow_text = self.tb_minx.get("1.0", "end-1c")
        xhigh_text = self.tb_maxx.get("1.0", "end-1c")
        
        try:
            xlow = float(xlow_text)
            xhigh = float(xhigh_text)
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter valid numeric values for X range.")
            return
        
        if xlow >= xhigh:
            messagebox.showerror("Error", "Minimum X value must be less than maximum X value.")
            return

        self.xrng = [xlow, xhigh]

# Set Y-Range Function
    def setYrng(self):
        ylow_text = self.tb_miny.get("1.0", "end-1c")
        yhigh_text = self.tb_maxy.get("1.0", "end-1c")
        
        try:
            ylow = float(ylow_text)
            yhigh = float(yhigh_text)
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter valid numeric values for Y range.")
            return
        
        if ylow >= yhigh:
            messagebox.showerror("Error", "Minimum Y value must be less than maximum Y value.")
            return

        self.yrng = [ylow, yhigh]


    def resetYrng(self,proj):
        # Only needed in undo, and only if what you want to
        # undo changed the y axis
        if ("setVelocity" in proj.history[-1]) or ("topoCorrect" in proj.history[-1]) and not self.picking: 
            self.yrng=self.prevyrng

    #Set Aspect Ratio Function
    def setAspect(self):
        aspect_ratio_text = self.tb_aspr.get("1.0", "end-1c")
        try:
            self.asp = float(aspect_ratio_text)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for aspect ratio. Please enter a valid number.")


    #######################################################################################
    #Velocity Controls functions
    # def setVelocity(self,proj):
    #     velocity =  sd.askfloat("Input","Radar wave velocity [m/ns]?")        
    #     if velocity is not None:
    #         proj.setVelocity(velocity)
    #         self.prevyrng=self.yrng
    #         self.yrng=[0,np.max(proj.depth)]

    def setVelocity(self,proj):
        velocity_str = self.tb_velo.get("1.0", "end-1c")
        try:
            velocity = float(velocity_str)
            # Assuming proj is defined elsewhere in your code
            proj.setVelocity(velocity)
            self.prevyrng = self.yrng
            self.yrng = [0, np.max(proj.depth)]
        except ValueError:
            tk.messagebox.showerror("Error", "Invalid input! Please enter a valid number.")


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
        if topofile != '':
            out = self.getDelimiter()    
            proj.topoCorrect(topofile,self.delimiter)
            self.prevyrng=self.yrng
            self.yrng=[proj.minTopo-np.max(proj.depth),proj.maxTopo]


    # def topoCorrect(self, proj):
    #     topofile = self.tb_topo.get("1.0", "end-1c")  # Get the value from the text box
    #     if topofile:
    #         if proj.velocity is None:
    #             mesbox.showinfo("Topo Correct Error", "You have to set the velocity first")
    #             return
    #         out = self.getDelimiter()
    #         proj.topoCorrect(topofile, self.delimiter)
    #         self.prevyrng = self.yrng
    #         self.yrng = [proj.minTopo - np.max(proj.depth), proj.maxTopo]
            
        

    #######################################################################################################################
    #Profile Controls functions

    # AdjProfile to take input from text boxes
    def adjProfile(self, proj):
        flipit = messagebox.askyesno("Question", "Flip the profile (left to right)?")
        if flipit:
            proj.flipProfile()

        # Get values from the text boxes
        minPos = self.tb_adjprfMin.get("1.0", tk.END).strip()
        maxPos = self.tb_adjprfMax.get("1.0", tk.END).strip()

        if not minPos or not maxPos:
            messagebox.showerror("Error", "Please fill in both minimum and maximum positions.")
            return

        try:
            minPos = float(minPos)
            maxPos = float(maxPos)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for position.")
            return

        proj.adjProfile(minPos=minPos, maxPos=maxPos)
        self.xrng = [minPos, maxPos]

    # Set Zero Time Function
    def setZeroTime(self, proj):
        # Get the value from the textbox
        newZeroTime_str = self.tb_zt.get("1.0", tk.END).strip()

        # Check if the input is valid and convert it to float
        try:
            newZeroTime = float(newZeroTime_str)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for new zero time.")
            return

        # Perform the action with the new zero time
        proj.setZeroTime(newZeroTime=newZeroTime)

    # Truncate Y Function
    def truncateY(self, proj):
        maxY = self.tb_adjprf.get("1.0", "end-1c").strip()  # Get input from text box
        if not maxY:
            messagebox.showerror("Error", "Please fill in the maximum Y value.")
            return

        try:
            maxY = float(maxY)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for maximum Y value.")
            return

        proj.truncateY(maxY)

    # Cut profile
    def cut(self, proj):
        minX = self.tb_cutx.get("1.0", "end-1c").strip()
        maxX = self.tb_cuty.get("1.0", "end-1c").strip()

        if not minX or not maxX:
            messagebox.showerror("Error", "Please fill in both minimum and maximum X values.")
            return

        try:
            minX = float(minX)
            maxX = float(maxX)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for X values.")
            return

        proj.cut(minX, maxX)

    # Dewow Function
    def dewow(self, proj):
        # Get the value from the text box
        window_value = self.tb_dewow.get("1.0", "end-1c").strip()

        # Check if the input value is a valid numeric value
        if not window_value.isdigit():
            messagebox.showerror("Error", "Please enter a valid numeric value for window width.")
            return

        window = float(window_value)

        proj.dewow(window=window)

    # Remove Mean Trace Function
    def remMeanTrace(self, proj):
        ntraces = self.tb_rmt.get("1.0", "end-1c").strip()
        if not ntraces:
            messagebox.showerror("Error", "Please fill in the number of traces.")
            return

        try:
            ntraces = float(ntraces)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for the number of traces.")
            return

        proj.remMeanTrace(ntraces=ntraces)

    # Smooth Function
    def smooth(self, proj):
        window = self.tb_smth1.get("1.0", "end-1c").strip()
        if not window:
            messagebox.showerror("Error", "Please fill in the window value.")
            return

        try:
            window = float(window)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for window value.")
            return

        proj.smooth(window=window)

    # Profile Smooth Function
    def profileSmooth(self, proj):
        ntraces = self.tb_trnum.get("1.0", "end-1c").strip()
        noversample = self.tb_cpynum.get("1.0", "end-1c").strip()

        if not ntraces or not noversample:
            messagebox.showerror("Error", "Please fill in both the number of traces and the oversample value.")
            return

        try:
            ntraces = int(ntraces)
            noversample = int(noversample)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for number of traces or oversample value.")
            return

        proj.profileSmooth(ntraces, noversample)

    # Tpow Gain Function
    def tpowGain(self, proj):
        power = self.tb_smth.get("1.0", "end-1c").strip()
        if not power:
            messagebox.showerror("Error", "Please fill in a value.")
            return

        try:
            power = float(power)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for power value.")
            return

        proj.tpowGain(power=power)

    # AGC Gain Function
    def agcGain(self, proj):
        window = self.tb_agc.get("1.0", "end-1c").strip()
        if not window:
            messagebox.showerror("Error", "Please provide a value.")
            return

        try:
            window = float(window)
        except ValueError:
            messagebox.showerror("Error", "Invalid input.")
            return

        proj.agcGain(window=window)

    # Show hyperbola
    def showHyp(self, proj, a):
        x0_text = self.tb_hypbcent.get("1.0", "end-1c").strip()  # Remove leading/trailing whitespace
        t0_text = self.tb_hypbapx.get("1.0", "end-1c").strip()
        v_text = self.tb_hypbvelo.get("1.0", "end-1c").strip()

        if not x0_text or not t0_text or not v_text:
            messagebox.showerror("Error", "Please fill in all the hyperbola parameters.")
            return

        try:
            x0 = float(x0_text)
            t0 = float(t0_text)
            v = float(v_text)
        except ValueError:
            messagebox.showerror("Error", "Invalid input for hyperbola parameters.")
            return

        y = proj.profilePos - x0
        d = v * t0 / 2.0
        k = np.sqrt(d ** 2 + np.power(y, 2))
        t2 = 2 * k / v
        a.plot(proj.profilePos, t2, '--c', linewidth=3)
        self.hypx = x0
        self.hypt = t0
        self.hypv = v


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
	
	root = ctk.CTk() #; root.iconbitmap()

	GPRPyApp(root)

	root.mainloop()
	
			  
	

if __name__ == '__main__':
	main()
			


###### If in future you want to use matplotlib and vtk together
    # def load_vtk_file(self):
    #     # Open a file dialog to select the VTK or VTS file
    #     self.vtk_file = filedialog.askopenfilename(filetypes=[("VTK files", "*.vtk"), ("VTS files", "*.vts")])
    #     if self.vtk_file:
    #         # Call function to read and display the VTK file
    #         self.display_vtk_file()

    # def display_vtk_file(self):
    #     vtk_dataset = self.read_vtk_file(self.vtk_file)
        
    #     # Convert VTK dataset to numpy array
    #     np_data = vtk_to_numpy(vtk_dataset.GetPointData().GetScalars())
    #     dims = vtk_dataset.GetDimensions()
    #     np_data = np_data.reshape(dims, order='F')  # Reshape according to VTK array layout

    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')

    #     # Basic visualization: plotting slices as 2D planes in 3D space
    #     for i in range(dims[2]):  # Iterate over Z-axis slices
    #         slice = np_data[:, :, i]
            
    #         x, y = np.mgrid[0:slice.shape[0], 0:slice.shape[1]]
    #         z = np.full(slice.shape, i)
            
    #         # Plotting the slice
    #         ax.contourf(x, y, z, zdir='z', levels=[i-0.5, i+0.5], cmap="viridis", alpha=0.5)
    #         ax.contourf(x, y, slice, zdir='z', offset=i, cmap="viridis")

    #     ax.set_xlabel('X axis')
    #     ax.set_ylabel('Y axis')
    #     ax.set_zlabel('Z axis')
    #     plt.show()
        
    # def read_vtk_file(self, filename):
    #     # Determines file extension to choose appropriate reader
    #     if filename.endswith(".vtk"):
    #         reader = vtk.vtkStructuredPointsReader()
    #     elif filename.endswith(".vts"):
    #         reader = vtk.vtkXMLStructuredGridReader()
    #     else:
    #         raise ValueError("Unsupported file format")
    #     reader.SetFileName(filename)
    #     reader.Update()
    #     return reader.GetOutput()

    # def extract_slice(self, vtk_dataset, slice_index):
    #     # Assumes the slice is along the Z-axis for simplicity
    #     dims = vtk_dataset.GetDimensions()
    #     if slice_index < 0 or slice_index >= dims[2]:
    #         raise ValueError("Slice index out of range.")

    #     # Extracting the slice (example for structured points, adjust if necessary)
    #     slice_data = vtk_to_numpy(vtk_dataset.GetPointData().GetScalars())
    #     slice_shape = (dims[1], dims[0])  # Assuming XY plane slices
    #     slice_array = slice_data[slice_index * dims[0] * dims[1]:(slice_index + 1) * dims[0] * dims[1]]
    #     slice_array = slice_array.reshape(slice_shape)

    #     return slice_array

    # def plot_slice(self, slice_array):
    #     plt.imshow(slice_array, cmap='viridis', origin='lower')
    #     plt.colorbar()
    #     plt.show()

    # def prompt_export_data_points(self):
    #     filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
    #     if filepath:
    #         self.export_data_points(filepath)