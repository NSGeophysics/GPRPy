import os
import pickle
from tkinter.filedialog import SaveAs
from turtle import fd
from matplotlib.backend_tools import ToolSetCursor
import numpy as np
import matplotlib.pyplot as plt
import copy
from PyQt5.QtWidgets import QFileDialog, QInputDialog, QMessageBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from gprpy.toolbox import gprIO_DT1, gprIO_DZT
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pyevtk.hl import gridToVTK


# Variables for storing GPR data and related information
data = None
info = None
profilePos = None
twtt = None
history = ["mygpr = gp.gprpyCW()"]
dtype = None
linStAmp = None
hypStAmp = None
lins = list()
hyps = list()

class FileOperations:
    def __init__(self,filename=None,dtype=None):
        '''
        Initialization for a gprpyCW object. Initialization can be 
        empty or with a provided filename for the GPR data and 
        a data type.

        INPUT:
        filename     data file name. Currently supported formats:
                     .gpr (GPRPy), .DT1 (SnS), .DZT (GSSI), .rd3 (MALA),
                     and ENVI standard BSQ.
        dtype        data type. Either "CMP" or "WARR"
        '''
        # Inheriting the initializer from the gprpyProfile class
        super().__init__()
        self.history = ["mygpr = gp.gprpyCW()"]
        # Initialize previous for undo
        self.previous = {}
        self.dtype = dtype
        # Stacked amplitude plots
        self.linStAmp = None
        self.hypStAmp = None
        # Picked lines and hyperbolae
        self.lins = list()
        self.hyps = list()
        
        if (filename is not None) and (dtype is not None):
            self.importdata(filename,dtype)


    def init_previous(self):
        '''
        Initialization of data structure that contains the step 
        before the most recent action.
        '''
        previous = {}
        previous["data"] = data
        previous["twtt"] = twtt 
        previous["info"] = info
        previous["profilePos"] = profilePos
        previous["velocity"] = None
        previous["depth"] = None
        previous["maxTopo"] = None
        previous["minTopo"] = None
        previous["threeD"] = None
        previous["data_pretopo"] = None
        previous["twtt_pretopo"] = None
        histsav = copy.copy(history)
        previous["history"] = histsav
        return previous

    def openFileDialog(self):
        try:
            # Open file dialog to select a file
            filename, _ = QFileDialog.getOpenFileName(None, "Select File", "", "All Files (*);;DT1 Files (*.dt1)")

            if filename:
                dtype = self.get_data_type_from_user()
                # Call importdata method with the selected filename
                if dtype:
                    # Call importdata method with the selected filename and data type
                    self.importdata(filename)
        except Exception as e:
            print("An error occurred:", e)

            
    def get_data_type_from_user(self):
        # Options for the data type
        options = ["CMP", "WARR"]
        # Prompt user for the data type
        dtype, ok = QInputDialog.getItem(None, "Select Data Type", "Select Data Type:", options, 0, False)
        if ok:
            return dtype
        else:
            # If user cancels, return None
            return None
        
    def importdata(self, filename):
        file_name, file_ext = os.path.splitext(filename)
        
        try:
            if file_ext == ".DT1" or file_ext == ".HD" or file_ext == ".dt1" or file_ext == ".hd":
                if file_ext == ".DT1" or  file_ext == ".HD":
                    self.data = gprIO_DT1.readdt1(file_name + ".DT1")
                    self.info = gprIO_DT1.readdt1Header(file_name + ".HD")  
                else:
                    self.data = gprIO_DT1.readdt1(file_name + ".dt1")
                    self.info = gprIO_DT1.readdt1Header(file_name + ".hd")
                
                self.profilePos = np.linspace(self.info["Start_pos"], self.info["Final_pos"], self.info["N_traces"])

                sec_per_samp = self.info["Total_time_window"] / self.info["N_pts_per_trace"]
                tshift = self.info["TZ_at_pt"] * sec_per_samp
                
                self.twtt = np.linspace(0, self.info["Total_time_window"], self.info["N_pts_per_trace"]) - tshift

                self.antsep = self.info["Antenna_sep"]
                self.velocity = None
                self.depth = None
                self.maxTopo = None
                self.minTopo = None
                self.threeD = None
                self.data_pretopo = None
                self.twtt_pretopo = None
                
                self.init_previous()
                
                histstr = "mygpr.importdata('%s')" % (filename)
                self.history.append(histstr)                                
                
            elif file_ext == ".DZT":
                self.data, self.info = gprIO_DZT.readdzt(filename)

                if self.info["rhf_spm"] != 0:
                    self.profilePos = self.info["rhf_position"] + np.linspace(0.0, self.data.shape[1] / self.info["rhf_spm"], self.data.shape[1])
                else:
                    self.profilePos = self.info["rhf_position"] + np.linspace(0.0, self.data.shape[1] / self.info["rhf_sps"], self.data.shape[1])
                    
                self.twtt = np.linspace(0, self.info["rhf_range"], self.info["rh_nsamp"])

                self.antsep = 0
                self.velocity = None
                self.depth = None
                self.maxTopo = None
                self.minTopo = None
                self.threeD = None
                self.data_pretopo = None
                self.twtt_pretopo = None
                
                self.init_previous()
                
                histstr = "mygpr.importdata('%s')" % (filename)
                self.history.append(histstr)
        
            elif file_ext == ".gpr":
                with open(filename, 'rb') as f:
                    data, info, profilePos, twtt, history, antsep, velocity, depth, maxTopo, minTopo, threeD, data_pretopo, twtt_pretopo = pickle.load(f)
                self.data = data
                self.info = info
                self.profilePos = profilePos
                self.twtt = twtt
                self.history = history
                self.antsep = antsep
                self.velocity = velocity
                self.depth = depth
                self.maxTopo = maxTopo
                self.minTopo = minTopo
                self.threeD = threeD
                self.data_pretopo = data_pretopo
                self.twtt_pretopo = twtt_pretopo
                
                self.init_previous()
                
            else:
                print("Can only read dt1, DT1, hd, HD, DZT, dat, GPRhdr, rad, rd3, rd7, and gpr files")
            
            # Show success message box
            
            QMessageBox.information(None, "Success", "File uploaded successfully.")
            self.displayDataInViewingWindow()
            
        except Exception as e:
            print("An error occurred:", e)
            QMessageBox.warning(None, "Error", f"An error occurred: {str(e)}")

    def save(self,filename):
        '''
        Saves the processed data together with the processing and visualization
        history. Warning: The history stored in this file will contain the full 
        path to the file.

        INPUT:
        filename       name for .gpr file
        '''
        # Saving the objects:
        # Want to force the file name .gpr
        file_name, file_ext = os.path.splitext(filename)
        if not(file_ext=='.gpr'):
            filename = filename + '.gpr'
        with open(filename, 'wb') as f:  
            pickle.dump([self.data, self.info, self.profilePos, self.twtt,
                         self.history, self.antsep, self.velocity, self.depth,
                         self.maxTopo, self.minTopo, self.threeD, self.data_pretopo,
                         self.twtt_pretopo], f)            
        print("Saved " + filename)
        # Add to history string
        histstr = "mygpr.save('%s')" %(filename)
        self.history.append(histstr)

    def showHistory(self):
        '''
        Prints out processing and visualization history of a data set. 
        '''
        for i in range(0,len(self.history)):
            print(self.history[i])

    def writeHistory(self,proj):        
        filename = self.asksaveasfilename(defaultextension=".py")
        if filename is not '':
            proj.writeHistory(filename)
            print("Wrote history to " + filename)
            
    def exportVTK(self,outfile,gpsinfo,delimiter=',',thickness=0,aspect=1.0,smooth=True, win_length=51, porder=3):
        '''
        Turn processed profile into a VTK file that can be imported in 
        Paraview or MayaVi or other VTK processing and visualization tools.

        If three-dimensional topo information is provided (X,Y,Z or 
        Easting, Northing, Elevation), then the profile will be exported 
        in its three-dimensional shape.

        INPUT:
        outfile       file name for the VTK file
        gpsinfo       EITHER: n x 3 matrix containing x, y, and z or 
                              Easting, Northing, Elevation information
                      OR: file name for ASCII text file containing this
                          information
        delimiter     if topo file is provided: delimiter (by comma, or by tab)
                      [default: ',']. To set tab: delimiter='\t' 
        thickness     If you want your profile to be exported as a 
                      three-dimensional band with thickness, enter thickness
                      in meters [default: 0]
        aspect        aspect ratio in case you want to exaggerate z-axis.
                      default = 1. I recommend leaving this at 1 and using 
                      your VTK visualization software to set the aspect for
                      the representation.
        smooth        Want to smooth the profile's three-dimensional alignment
                      instead of piecewise linear? [Default: True]
        win_length    If smoothing, the window length for 
                      scipy.signal.savgol_filter [default: 51]
        porder        If smoothing, the polynomial order for
                      scipy.signal.savgol_filter [default: 3]
        '''
        # If gpsmat is a filename, we first need to load the file:
        if type(gpsinfo) is str:
            gpsmat = np.loadtxt(gpsinfo,delimiter=delimiter)
        else:
            gpsmat = gpsinfo
            
        # First get the x,y,z positions of our data points
        x,y,z = ToolSetCursor.prepVTK(self.profilePos,gpsmat,smooth,win_length,porder)        
        z = z*aspect     
        if self.velocity is None:
            downward = self.twtt*aspect
        else:
            downward = self.depth*aspect                        
        Z = np.reshape(z,(len(z),1)) - np.reshape(downward,(1,len(downward)))

        
        if thickness:
            ZZ = np.tile(np.reshape(Z, (1,Z.shape[0],Z.shape[1])), (2,1,1))
        else:
            ZZ = np.tile(np.reshape(Z, (1,Z.shape[0],Z.shape[1])), (1,1,1))
        
        # This is if we want everything on the x axis.
        #X = np.tile(np.reshape(self.profilePos,(len(self.profilePos),1)),(1,len(downward)))
        #XX = np.tile(np.reshape(X, (X.shape[0],1,X.shape[1])), (1,2,1))
        #YY = np.tile(np.reshape([-thickness/2,thickness/2],(1,2,1)), (len(x),1,len(downward)))

        # To create a 3D grid with a width, calculate the perpendicular direction,
        # normalize it, and add it to xvals and yvals as below.
        # To figure this out, just drar the profile point-by-point, and at each point,
        # draw the perpendicular to the segment and place a grid point in each perpendicular
        # direction
        #
        #          x[0]-px[0], x[1]-px[1], x[2]-px[2], ..... 
        # xvals =     x[0]   ,    x[1]   ,     x[2]  , .....   
        #          x[0]+px[0], x[1]+px[1], x[2]+px[2], .....
        #  
        #          y[0]+py[0], y[1]+py[1], y[2]+py[2], .....
        # yvals =     y[0]   ,    y[1]   ,    y[2]   , .....
        #          y[0]-py[0], y[1]-py[1], y[2]-py[2], .....
        #
        # Here, the [px[i],py[i]] vector needs to be normalized by the thickness
        if thickness:
            pvec = np.asarray([(y[0:-1]-y[1:]).squeeze(), (x[1:]-x[0:-1]).squeeze()])
            pvec = np.divide(pvec, np.linalg.norm(pvec,axis=0)) * thickness/2.0
            # We can't calculate the perpendicular direction at the last point
            # let's just set it to the same as for the second-to-last point
            pvec = np.append(pvec, np.expand_dims(pvec[:,-1],axis=1) ,axis=1)
            X = np.asarray([(x.squeeze()-pvec[0,:]).squeeze(), (x.squeeze()+pvec[0,:]).squeeze()])
            Y = np.asarray([(y.squeeze()+pvec[1,:]).squeeze(), (y.squeeze()-pvec[1,:]).squeeze()])
        else:
            X = np.asarray([x.squeeze()])
            Y = np.asarray([y.squeeze()])
        
        # Copy-paste the same X and Y positions for each depth
        XX = np.tile(np.reshape(X, (X.shape[0],X.shape[1],1)), (1,1,ZZ.shape[2]))
        YY = np.tile(np.reshape(Y, (Y.shape[0],Y.shape[1],1)), (1,1,ZZ.shape[2]))
        
        if self.maxTopo is None:
            data=self.data.transpose()
        else:
            data=self.data_pretopo.transpose()       

        data = np.asarray(data)
        data = np.reshape(data,(1,data.shape[0],data.shape[1]))                 
        data = np.tile(data, (2,1,1))
        
        # Remove the last row and column to turn it into a cell
        # instead of point values 
        data = data[0:-1,0:-1,0:-1]

        nx=2-1
        ny=len(x)-1
        nz=len(downward)-1
        datarray = np.zeros(nx*ny*nz).reshape(nx,ny,nz)
        datarray[:,:,:] = data
        
        gridToVTK(outfile, XX, YY, ZZ, cellData={'gpr': datarray})
 
        # Put in history
        if gpsinfo is None:
            histstr = "mygpr.exportVTK('%s',aspect=%g)" %(outfile,aspect)
        else:
            if type(gpsinfo) is str:            
                if delimiter is ',':
                    histstr = "mygpr.exportVTK('%s',gpsinfo='%s',thickness=%g,delimiter=',',aspect=%g,smooth=%r, win_length=%d, porder=%d)" %(outfile,gpsinfo,thickness,aspect,smooth,win_length,porder)
                else:
                    histstr = "mygpr.exportVTK('%s',gpsinfo='%s',thickness=%g,delimiter='\\t',aspect=%g,smooth=%r, win_length=%d, porder=%d)" %(outfile,gpsinfo,thickness,aspect,smooth,win_length,porder)
            else:
                if delimiter is ',':
                    histstr = "mygpr.exportVTK('%s',gpsinfo=mygpr.threeD,thickness=%g,delimiter=',',aspect=%g,smooth=%r, win_length=%d, porder=%d)" %(outfile,thickness,aspect,smooth,win_length,porder)
                else:
                    histstr = "mygpr.exportVTK('%s',gpsinfo=mygpr.threeD,thickness=%g,delimiter='\\t',aspect=%g,smooth=%r, win_length=%d, porder=%d)" %(outfile,thickness,aspect,smooth,win_length,porder)
                    
        self.history.append(histstr)   



    def prepCWFig(self, contrast=1.0, color="gray", yrng=None, xrng=None, showlnhp=False):
        '''
        This is a helper function.
        It prepares the plot showing the processed CMP or WARR data.
        
        INPUT:
        color        "gray", or "bwr" for blue-white-red,
                        or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                        [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        showlnhp     show the observed lines and hyperbolae from the list
                        [default: False]

        OUTPUT:
        contrast     contrast value used to prepare the figure 
        color        color value used to prepare the figure
        yrng         yrng value used to prepare the figure
        xrng         xrng value used to prepare the figure 
        showlnhp     showlnhp value used to prepare the figure
        '''
        dx=self.profilePos[3]-self.profilePos[2]
        dt=self.twtt[3]-self.twtt[2]
        stdcont = np.nanmax(np.abs(self.data)[:])       
        
        plt.imshow(self.data,cmap=color,extent=[min(self.profilePos)-dx/2.0,
                                                max(self.profilePos)+dx/2.0,
                                                max(self.twtt)+dt/2.0,
                                                min(self.twtt)-dt/2.0],
                    aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)
        plt.gca().set_ylabel("time [ns]")
        plt.gca().invert_yaxis()
        if yrng is not None:
            yrng=[np.max(yrng),np.min(yrng)]
        else:
            yrng=[np.max(self.twtt),np.min(self.twtt)]
        plt.ylim(yrng)
            
        if xrng is None:
            xrng=[min(self.profilePos),max(self.profilePos)]                           
        plt.xlim(xrng)

        plt.gca().get_xaxis().set_visible(True)
        plt.gca().get_yaxis().set_visible(True)
        if self.dtype == "WARR":
            plt.gca().set_xlabel("antenna separation [m]")
            typefact=1
        elif self.dtype == "CMP":
            plt.gca().set_xlabel("distance from midpoint [m]")
            typefact=2
        plt.gca().xaxis.tick_top()
        plt.gca().xaxis.set_label_position('top')

        # Show hyperbolae / lines if you want
        if showlnhp:
            if self.lins:
                for lin in self.lins:
                    time = lin[0] + typefact*self.profilePos/lin[1]
                    plt.plot(self.profilePos,time,linewidth=2,color='yellow')
                    plt.plot(self.profilePos,time,linewidth=1,color='black')
            if self.hyps:
                x2 = np.power(typefact*self.profilePos,2.0)
                for hyp in self.hyps:
                    time = np.sqrt(x2 + 4*np.power(hyp[0]/2.0 * hyp[1],2.0))/hyp[1]
                    plt.plot(self.profilePos,time,linewidth=2,color='yellow')
                    plt.plot(self.profilePos,time,linewidth=1,color='black')
                                    
        return contrast, color, yrng, xrng, showlnhp


    def prepStAmpFig(self, whichstamp="lin", saturation=1.0, yrng=None, vrng=None):
        '''
        This is a helper function.
        It prepares the plot showing the stacked amplitudes results.
        
        INPUT:
        whichstamp   is this for the linear ("lin") or hyperbolic ("hyp") 
                        stacked amplitudes
        saturation   Factor to increase contrast by reducing color range.
                        [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        vrng         velocities (x-axis) range to show 
                        [default: None, meaning "everything"]

        OUTPUT:
        whichstamp   whichstamp value used to prepare the figure
        saturation   saturation value used to prepare the figure 
        yrng         yrng value used to prepare the figure 
        vrng         vrng value used to prepare the figure
        '''
        dt=self.twtt[3]-self.twtt[2]
        dv=self.vVals[3]-self.vVals[2]
        if whichstamp == "lin":
            stamp = self.linStAmp
            title = "linear stacked amplitude"
        elif whichstamp == "hyp":
            stamp = self.hypStAmp
            title = "hyperbolic stacked amplitude"
        else:
            stamp = None
            
        if stamp is not None:
            stdcont = np.nanmax(np.abs(stamp)[:])
            plt.imshow(np.flipud(np.abs(stamp)), cmap='inferno',
                        extent=[np.min(self.vVals)-dv/2.0,
                                np.max(self.vVals)+dv/2.0,
                                np.min(self.twtt)-dt/2.0,
                                np.max(self.twtt)+dt/2.0],
                        aspect='auto',
                        vmin=0, vmax=stdcont/saturation)

            if yrng is not None:
                yrng=[np.max(yrng),np.min(yrng)]
            else:
                yrng=[np.max(self.twtt),np.min(self.twtt)]
            plt.ylim(yrng)
            
            if vrng is None:
                vrng=[np.min(self.vVals),np.max(self.vVals)]                           
            plt.xlim(vrng)
                

            plt.gca().set_xlabel("velocity [m/ns]")
            plt.gca().set_ylabel("time [ns]")
            #plt.gca().invert_yaxis()
            plt.gca().set_title(title)
            plt.gca().get_xaxis().set_visible(True)
            plt.gca().get_yaxis().set_visible(True)
            plt.gca().get_xaxis().set_ticks_position('both')
            plt.gca().get_yaxis().set_ticks_position('both')
                                
        return whichstamp, saturation, yrng, vrng

    def showCWFig(self, **kwargs):
        '''
        Plots the CMP or WARR data using Matplotlib. 
        You need to run .show() afterward to show it 

        INPUT:
        color        "gray", or "bwr" for blue-white-red,
                        or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                        [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        showlnhp     show the observed lines and hyperbolae from the list
                        [default: False]
        '''
        self.prepCWFig(**kwargs)
        plt.show(block=False)


    def showStAmpFig(self, **kwargs):
        '''
        Plots the stacked amplitude results using Matplotlib. 
        You need to run .show() afterward to show it 

        INPUT:
        whichstamp   is this for the linear ("lin") or hyperbolic ("hyp") 
                        stacked amplitudes
        saturation   Factor to increase contrast by reducing color range.
                        [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        vrng         velocities (x-axis) range to show 
                        [default: None, meaning "everything"]
        '''
        self.prepStAmpFig(**kwargs)
        plt.show(block=False)


    def displayDataInViewingWindow(self):
        # Create a canvas to display the plots
        canvas = FigureCanvas(plt.figure())

        # Call the prepCWFig method to prepare the plot
        self.prepCWFig()

        # Plot the CMP or WARR data on the canvas
        self.showCWFig()

        # If you want to plot stacked amplitude results, uncomment the following lines:
        # self.prepStAmpFig()
        # self.showStAmpFig()

        # Return the canvas object
        return canvas

    def asksaveasfilename(**options):
        "Ask for a filename to save as"
        return SaveAs(**options).show()



