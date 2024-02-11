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

class GPRPyApp:
    def __init__(self, master):

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

        self.window = master
        master.geometry('1024x768')
        master.title("GPR - Py")
        master.rowconfigure(0, minsize=766, weight=1)
        master.columnconfigure(1, minsize=766, weight=1)
       

        proj = gp.gprpyProfile()

        btn_frm = tk.Frame(master, relief= tk.RAISED, bd = 2)#, height= 700, width= 100)
        btn_frm.grid(row = 0, column= 0, sticky= 'ns')


        fig=Figure(figsize=(self.widfac,self.highfac))
        a=fig.add_subplot(111)
        mpl.rcParams.update({'font.size': mpl.rcParams['font.size']*self.widfac})
        a.tick_params(direction='out',length=6*self.widfac,width=self.highfac)
        
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.get_tk_widget().grid(row=2,column=0,columnspan= 9,rowspan= 22,sticky='nsew')

        canvas.draw() 

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

# Button and checkbutton, these will
# appear in collapsible pane container

def main():
	rightcol=9
	figrowsp=19+1
	
	root = tk.Tk()
	GPRPyApp(root)
	root.mainloop()
	
			  
	

if __name__ == '__main__':
	main()
			