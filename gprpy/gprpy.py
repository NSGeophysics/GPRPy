import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gprpy.toolbox.gprIO_DT1 as gprIO_DT1
import gprpy.toolbox.gprIO_DZT as gprIO_DZT
import gprpy.toolbox.gprIO_BSQ as gprIO_BSQ
import gprpy.toolbox.gprIO_MALA as gprIO_MALA
import gprpy.toolbox.gprpyTools as tools
try:
    import gprpy.irlib.external.mig_fk as mig_fk
except:
    print("Install fk migration if needed")
import copy
import scipy.interpolate as interp
from pyevtk.hl import gridToVTK

class gprpyProfile:
    '''
    Ground penetrating radar data processing and visualization class 
    for common-offset profiles.
    '''

    def __init__(self,filename=None):
        '''
        Initialization for a gprpyProfile object. Initialization can be 
        empty or with a provided filename for the GPR data.

        INPUT:
        filename     data file name. Currently supported formats:
                     .gpr (GPRPy), .DT1 (SnS), .DZT (GSSI), .rd3 (MALA),
                     and ENVI standard BSQ.
        '''
        
        self.history = ["mygpr = gp.gprpyProfile()"]

        # Initialize previous for undo
        self.previous = {}
        
        if filename is not None:
            self.importdata(filename)                 
        
    def importdata(self,filename):
        '''
        Loads .gpr (native GPRPy), .DT1 (Sensors and Software),
        .DZT (GSSI), .GPRhdr (ENVI standard BSQ), .rad (MALA)
        data files and populates all the gprpyProfile fields.

        INPUT: 
        filename  name of the .gpr, .DT1, dt1, .DZT, .GPRhdr, dat, 
                  rd3, or .rad file you want to import.
                  The header file name and the data file name 
                  have to be the same!
        '''
        
        file_name, file_ext = os.path.splitext(filename)
        
        if file_ext==".DT1" or file_ext==".HD" or file_ext==".dt1" or file_ext==".hd":
            if file_ext==".DT1" or  file_ext==".HD":
                self.data=gprIO_DT1.readdt1(file_name + ".DT1")
                self.info=gprIO_DT1.readdt1Header(file_name + ".HD")  
            else:
                self.data=gprIO_DT1.readdt1(file_name + ".dt1")
                self.info=gprIO_DT1.readdt1Header(file_name + ".hd")
            
            self.profilePos = np.linspace(self.info["Start_pos"],
                                          self.info["Final_pos"],
                                          self.info["N_traces"])

            #self.twtt = np.linspace(self.info["TZ_at_pt"],
            #                        self.info["Total_time_window"],
            #                        self.info["N_pts_per_trace"])

            sec_per_samp = self.info["Total_time_window"]/self.info["N_pts_per_trace"]
            tshift = self.info["TZ_at_pt"]*sec_per_samp
            
            self.twtt = np.linspace(0,self.info["Total_time_window"],
                                    self.info["N_pts_per_trace"]) - tshift

            self.antsep = self.info["Antenna_sep"] # Set to m in the loading routine 
            self.velocity = None
            self.depth = None
            self.maxTopo = None
            self.minTopo = None
            self.threeD = None
            self.data_pretopo = None
            self.twtt_pretopo = None
            # Initialize previous
            self.initPrevious()
            
            # Put what you did in history
            histstr = "mygpr.importdata('%s')" %(filename)
            self.history.append(histstr)                                
            
        elif file_ext==".DZT":

            self.data, self.info = gprIO_DZT.readdzt(filename)

            if self.info["rhf_spm"] != 0:
                self.profilePos = self.info["rhf_position"]+np.linspace(0.0,
                                                                        self.data.shape[1]/self.info["rhf_spm"],
                                                                        self.data.shape[1])
            else:
                self.profilePos = self.info["rhf_position"]+np.linspace(0.0,
                                                                        self.data.shape[1]/self.info["rhf_sps"],
                                                                        self.data.shape[1])
                
            self.twtt = np.linspace(0,self.info["rhf_range"],self.info["rh_nsamp"])

            self.antsep = 0
            self.velocity = None
            self.depth = None
            self.maxTopo = None
            self.minTopo = None
            self.threeD = None
            self.data_pretopo = None
            self.twtt_pretopo = None
            # Initialize previous
            self.initPrevious()
            
            # Put what you did in history
            histstr = "mygpr.importdata('%s')" %(filename)
            self.history.append(histstr)
    
        

        elif file_ext==".GPRhdr" or file_ext==".dat":
            # ENVI standard BSQ file
            self.data, self.info = gprIO_BSQ.readBSQ(file_name)

            self.profilePos = float(self.info["dx"])*np.arange(0,int(self.info["columns"]))
            self.twtt = np.linspace(0,float(self.info["time_window"]),int(self.info["lines"]))

            self.antsep = 0
            self.velocity = None
            self.depth = None
            self.maxTopo = None
            self.minTopo = None
            self.threeD = None
            self.data_pretopo = None
            self.twtt_pretopo = None
            # Initialize previous
            self.initPrevious()
            
            # Put what you did in history
            histstr = "mygpr.importdata('%s')" %(filename)
            self.history.append(histstr)       


        elif file_ext==".rad" or file_ext==".rd3" or file_ext==".rd7":
            self.data, self.info = gprIO_MALA.readMALA(file_name)

            self.twtt = np.linspace(0,float(self.info["TIMEWINDOW"]),int(self.info["SAMPLES"]))
            self.profilePos = float(self.info["DISTANCE INTERVAL"])*np.arange(0,self.data.shape[1])

            self.antsep = self.info["ANTENNA SEPARATION"]
            self.velocity = None
            self.depth = None
            self.maxTopo = None
            self.minTopo = None
            self.threeD = None
            self.data_pretopo = None
            self.twtt_pretopo = None
            # Initialize previous
            self.initPrevious()
            
            # Put what you did in history
            histstr = "mygpr.importdata('%s')" %(filename)
            self.history.append(histstr)
            
            

        elif file_ext==".gpr":
            ## Getting back the objects:
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
            
            # Initialize previous
            self.initPrevious()
            
        else:
            print("Can only read dt1, DT1, hd, HD, DZT, dat, GPRhdr, rad, rd3, rd7, and gpr files")

    def showHistory(self):
        '''
        Prints out processing and visualization history of a data set. 
        '''
        for i in range(0,len(self.history)):
            print(self.history[i])

    def writeHistory(self,outfilename="myhistory.py"):
        '''
        Turns the processing and visualization history into a Python script.
        The full path names are saved in the Python script. You can edit the
        Python script after saving to remove the full path names.

        INPUT:
        outfilename        filename for Python script
        '''
        with open(outfilename,"w") as outfile:
            outfile.write("# Automatically generated by GPRPy\nimport gprpy.gprpy as gp\n")
            for i in range(0,len(self.history)):
                outfile.write(self.history[i] + "\n")
                
    def undo(self):
        '''
        Undoes the last processing step and removes that step fromt he history.
        '''
        self.data = self.previous["data"]
        self.twtt = self.previous["twtt"]
        self.info = self.previous["info"]
        self.profilePos = self.previous["profilePos"]
        self.velocity = self.previous["velocity"]
        self.depth = self.previous["depth"]
        self.maxTopo = self.previous["maxTopo"]
        self.minTopo = self.previous["minTopo"]
        self.threeD = self.previous["threeD"]
        self.data_pretopo = self.previous["data_pretopo"]
        self.twtt_pretopo = self.previous["twtt_pretopo"]
        # Make sure to not keep deleting history
        # when applying undo several times. 
        histsav = copy.copy(self.previous["history"])
        del histsav[-1]
        self.history = histsav
        print("undo")

        
    def initPrevious(self):
        '''
        Initialization of data strucure that contains the step 
        before the most recent action.
        '''
        self.previous["data"] = self.data
        self.previous["twtt"] = self.twtt 
        self.previous["info"] = self.info
        self.previous["profilePos"] = self.profilePos
        self.previous["velocity"] = self.velocity
        self.previous["depth"] = self.depth
        self.previous["maxTopo"] = self.maxTopo
        self.previous["minTopo"] = self.minTopo
        self.previous["threeD"] = self.threeD
        self.previous["data_pretopo"] = self.data_pretopo
        self.previous["twtt_pretopo"] = self.twtt_pretopo
        histsav = copy.copy(self.history)
        self.previous["history"] = histsav

        

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

    
    # This is a helper function
    def prepProfileFig(self, color="gray", contrast=1.0, yrng=None, xrng=None, asp=None):
        '''
        This is a helper function.
        It prepares the plot showing the processed profile data.
        
        INPUT:
        color        "gray", or "bwr" for blue-white-red,
                     or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                     [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        asp          aspect ratio [default: None, meaning automatic]

        OUTPUT:
        contrast     contrast value used to prepare the figure
        color        color value used to prepare the figure
        yrng         yrng value used to prepare the figure
        xrng         xrng value used to prepare the figure
        asp          asp value used to prepare the figure
        ax           figure axis

        '''
        dx=self.profilePos[3]-self.profilePos[2]
        dt=self.twtt[3]-self.twtt[2]
        stdcont = np.nanmax(np.abs(self.data)[:])       
        
        if self.velocity is None:
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
            
        elif self.maxTopo is None:
            dy=dt*self.velocity
            plt.imshow(self.data,cmap=color,extent=[min(self.profilePos)-dx/2.0,
                                                    max(self.profilePos)+dx/2.0,
                                                    max(self.depth)+dy/2.0,
                                                    min(self.depth)-dy/2.0],
                       aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)
            plt.gca().set_ylabel("depth [m]")
            plt.gca().invert_yaxis()
            if yrng is not None:
                yrng=[np.max(yrng),np.min(yrng)]
            else:
                yrng=[np.max(self.depth),np.min(self.depth)]
                
        else:
            dy=dt*self.velocity
            plt.imshow(self.data,cmap=color,extent=[min(self.profilePos)-dx/2.0,
                                                    max(self.profilePos)+dx/2.0,
                                                    self.minTopo-max(self.depth)-dy/2.0,
                                                    self.maxTopo-min(self.depth)+dy/2.0],
                    aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)            
            plt.gca().set_ylabel("elevation [m]")
            if yrng is None:
                yrng=[self.minTopo-np.max(self.depth),self.maxTopo-np.min(self.depth)]
            

        if xrng is None:
            xrng=[min(self.profilePos),max(self.profilePos)]       
                
        if yrng is not None:
            plt.ylim(yrng)
            
        if xrng is not None:
            plt.xlim(xrng)

        if asp is not None:
            plt.gca().set_aspect(asp)

        plt.gca().get_xaxis().set_visible(True)
        plt.gca().get_yaxis().set_visible(True)                
        plt.gca().set_xlabel("profile position [m]")
        plt.gca().xaxis.tick_top()
        plt.gca().xaxis.set_label_position('top')
        ax = plt.gca()
        
        return contrast, color, yrng, xrng, asp, ax
       
    
    def showProfile(self, **kwargs):
        '''
        Plots the profile using Matplotlib. 
        You need to run .show() afterward to show it 
        
        INPUT:
        color        "gray", or "bwr" for blue-white-red,
                     or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                     [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        asp          aspect ratio [default: None, meaning automatic]

        '''
        contrast, color, yrng, xrng, asp, ax = self.prepProfileFig(**kwargs)
        #plt.show(block=False)
        return ax


    def printProfile(self, figname, dpi=600, **kwargs):
        '''
        Creates a pdf of the profile. 
        
        INPUT:
        figname      file name for the pdf
        dpi          dots per inch resolution [default: 600 dpi]
        color        "gray", or "bwr" for blue-white-red,
                     or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                     [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        asp          aspect ratio [default: None, meaning automatic]

        '''
        contrast, color, yrng, xrng, asp, ax = self.prepProfileFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        # Put what you did in history
        if asp is None:
            histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],dpi)
        else:
            histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], asp=%g, dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],asp,dpi)
        self.history.append(histstr)
        

    ####### Processing #######
    
    def adjProfile(self,minPos,maxPos):
        '''
        Adjusts the length of the profile.

        INPUT:
        minPos      starting position of the profile
        maxpos      end position of the profile
        '''
        # Store previous state for undo
        self.storePrevious()
        # set new profile positions
        self.profilePos = np.linspace(minPos,maxPos,len(self.profilePos))       
        # Put what you did in history
        histstr = "mygpr.adjProfile(%g,%g)" %(minPos,maxPos)
        self.history.append(histstr)

    
    def flipProfile(self):
        '''
        Flips the profile left-to-right (start to end)
        '''
        # Flips the profile left to right (start to end)
        self.storePrevious()
        self.data=np.flip(self.data,1)
        if self.data_pretopo is not None:
            self.data_pretopo = np.flip(self.data_pretopo,1)
        histstr = "mygpr.flipProfile()"
        self.history.append(histstr)
        

    def alignTraces(self):
        '''
        Aligns the traces in the profile such that their maximum 
        amplitudes align at the average travel time of the 
        maximum amplitudes.
        '''
        # Store previous state for undo
        self.storePrevious()        
        self.data = tools.alignTraces(self.data)      
        # Put what you did in history
        histstr = "mygpr.alignTraces()"
        self.history.append(histstr)



    def cut(self,minPos,maxPos):
        '''
        Removes all data outside of the profile positions between
        minPos and maxPos.

        INPUT:
        minPos      starting position of data to keep
        maxPos      end position of data to keep
        '''
        # Store previous state for undo
        self.storePrevious()
        zeroind = np.abs(self.profilePos - minPos).argmin()
        maxind = np.abs(self.profilePos - maxPos).argmin()
        self.data = self.data[:,zeroind:(maxind+1)]
        self.profilePos=self.profilePos[zeroind:(maxind+1)]
        if self.data_pretopo is not None:
            self.data_pretopo = self.data_pretopo[:,zeroind:(maxind+1)]
        # Put into history string
        histstr = "mygpr.cut(%g,%g)" %(minPos,maxPos)
        self.history.append(histstr)
        
        
    def setZeroTime(self,newZeroTime):
        '''
        Deletes all data recorded before newZeroTime and 
        sets newZeroTime to zero.

        INPUT:
        newZeroTime     The new zero-time
        '''
        # Store previous state for undo
        self.storePrevious()
        # Find index of value that is nearest to newZeroTime
        zeroind = np.abs(self.twtt - newZeroTime).argmin() 
        # Cut out everything before
        self.twtt = self.twtt[zeroind:] - newZeroTime
        # Set first value to 0
        self.twtt[0] = 0
        self.data = self.data[zeroind:,:]
        # Put what you did in history
        histstr = "mygpr.setZeroTime(%g)" %(newZeroTime)
        self.history.append(histstr)  

        
    def dewow(self,window):
        '''
        Subtracts from each sample along each trace an 
        along-time moving average.

        Can be used as a low-cut filter.

        INPUT:
        window     length of moving average window 
                   [in "number of samples"]
        '''
        # Store previous state for undo
        self.storePrevious()
        self.data = tools.dewow(self.data,window)
        # Put in history
        histstr = "mygpr.dewow(%d)" %(window)
        self.history.append(histstr)


    def smooth(self,window):
        '''
        Replaces each sample along each trace with an 
        along-time moving average.

        Can be used as high-cut filter.

        INPUT: 
        window     length of moving average window
                   [in "number of samples"]
        '''
        # Store previous state for undo
        self.storePrevious()
        self.data = tools.smooth(self.data,window)
        # Put in history
        histstr = "mygpr.smooth(%d)" %(window)
        self.history.append(histstr)

        
    def remMeanTrace(self,ntraces):
        '''
        Subtracts from each trace the average trace over
        a moving average window.

        Can be used to remove horizontal arrivals, 
        such as the airwave.

        INPUT:
        ntraces     window width; over how many traces 
                    to take the moving average. 
        '''
        # Store previous state for undo
        self.storePrevious()
        # apply
        self.data = tools.remMeanTrace(self.data,ntraces)        
        # Put in history
        histstr = "mygpr.remMeanTrace(%d)" %(ntraces)
        self.history.append(histstr)


    def profileSmooth(self,ntraces,noversample):
        '''
        First creates copies of each trace and appends the copies 
        next to each trace, then replaces each trace with the 
        average trace over a moving average window.

        Can be used to smooth-out noisy reflectors appearing 
        in neighboring traces, or simply to increase the along-profile 
        resolution by interpolating between the traces.

        For example: To increase the along-profile resolution smoothly 
        by a factor of 4: use

        mygpr.profileSmooth(4,4)

        INPUT:
        ntraces         window width [in "number of samples"]; 
                        over how many traces to take the moving average. 
        noversample     how many copies of each trace
        '''
        # Store previous state for undo
        self.storePrevious()
        self.data,self.profilePos = tools.profileSmooth(self.data,self.profilePos,
                                                        ntraces,noversample)
        # Put in history
        histstr = "mygpr.profileSmooth(%d,%d)" %(ntraces,noversample)
        self.history.append(histstr)

        
    def tpowGain(self,power=0.0):
        '''
        Apply a t-power gain to each trace with the given exponent.

        INPUT:
        power     exponent
        '''
        # Store previous state for undo
        self.storePrevious()
        # apply tpowGain
        self.data = tools.tpowGain(self.data,self.twtt,power)
        # Put in history
        histstr = "mygpr.tpowGain(%g)" %(power)
        self.history.append(histstr)

    def agcGain(self,window=10):
        '''
        Apply automated gain controll (AGC) by normalizing the energy
        of the signal over a given window width in each trace

        INPUT:
        window     window width [in "number of samples"]
        '''
        # Store previous state for undo
        self.storePrevious()
        # apply agcGain
        self.data = tools.agcGain(self.data,window)
        # Put in history
        histstr = "mygpr.agcGain(%d)" %(float(window))
        self.history.append(histstr)
        

    def setVelocity(self,velocity):
        '''
        Provide the subsurface RMS radar velocity

        INPUT:
        velocity      subsurface RMS velocity [in m/ns]
        '''
        # Store previous state for undo
        self.storePrevious()

        self.velocity = velocity
        self.depth = self.twtt * velocity/2.0
        
        # Put in history
        histstr = "mygpr.setVelocity(%g)" %(velocity)
        self.history.append(histstr)


    def antennaSep(self):
        ''' 
        Corrects for distortions of arrival times caused by the
        separation of the antennae.

        For this to work properly, you must have set the velocity
        and you must have set the zero time to the beginning of the 
        arrival of the airwave.

        '''

        # Store previous state for undo
        self.storePrevious()

        # Take into account that the airwave first break
        # is after the airwave has already traveled the
        # antenna separation with the speed of light 0.3 m/ns.
        # And we only look at half the
        # two-way travel time. Hence divide by two
        t0 = self.twtt/2 + self.antsep/(2*0.3)

        # t0 is when the waves left the transmitter antenna.
        # To be able to calculate the depth time from the
        # single-way travel time we need to shift the time reference
        # frame. Lets set it "arriving at midpoint time", so
        ta = t0 + self.antsep/(2*self.velocity)
        # Later we will need to undo this reference frame transformation

        # Now use the pythagorean relationship between single-way travel
        # time to the depth point and the depth time td
        tad = np.sqrt( ta**2 - (self.antsep/(2*self.velocity))**2 )

        # We are still in the "arriving at midpoint" time frame ta
        # To transform ta into depth time td, we need to shift it back
        # by the time it took for the ground wave to get to the midpoint.
        # This makes sense because the times before the groundwave got to the
        # midpoint will not actually be underground in the sense of:
        # No travel into depth has been recorded at the receiver.
        # These "arrivals" will just be shifted into "negative arrival times"
        # and hence "negative depth"
        td = tad - self.antsep/(2*self.velocity)

        # Finally, translate time into depth
        self.depth = td*self.velocity

        # And update the two-way travel time
        self.twtt = td
        
        # Put in history
        histstr = "mygpr.antennaSep()"
        self.history.append(histstr)

        
    def fkMigration(self):
        '''
        Apply Stolt's f-k migration to the profile. Requires the 
        velocity to be set.

        This is a wrapper function for the migration code
        imported from Nat Wilson's irlib software.
        '''
        # Store previous state for undo
        self.storePrevious()
        # apply migration
        dt=self.twtt[3]-self.twtt[2]
        #dx=self.profilePos[1]-self.profilePos[0]
        dx=(self.profilePos[-1]-self.profilePos[0])/(len(self.profilePos)-1)
        # fkmig sets x profile to start at zero but resamples
        self.data,self.twtt,migProfilePos=mig_fk.fkmig(self.data,dt,dx,self.velocity)
        self.profilePos = migProfilePos + self.profilePos[0]
        
        # Put in history
        histstr = "mygpr.fkMigration()"
        self.history.append(histstr)
        
        
    def truncateY(self,maxY):
        '''
        Delete all data after y-axis position maxY.

        INPUT:
        maxY    maximum y-axis position for data to be kept
        '''
        # Store previous state for undo
        self.storePrevious()
        if self.velocity is None:
            maxtwtt = maxY
            maxind = np.argmin( np.abs(self.twtt-maxY) )
            self.twtt = self.twtt[0:maxind]
            # Set the last value to maxY
            self.twtt[-1] = maxY
            self.data = self.data[0:maxind,:]
        else:
            maxtwtt = maxY*2.0/self.velocity
            maxind = np.argmin( np.abs(self.twtt-maxtwtt) )
            self.twtt = self.twtt[0:maxind]
            # Set the last value to maxtwtt
            self.twtt[-1] = maxtwtt
            self.data = self.data[0:maxind,:]
            self.depth = self.depth[0:maxind]
            self.depth[-1] = maxY
        # Put in history
        histstr = "mygpr.truncateY(%g)" %(maxY)
        self.history.append(histstr)


        
    def topoCorrect(self,topofile,delimiter=','):
        '''
        Correct for topography along the profile by shifting each 
        Trace up or down depending on a provided ASCII text file
        containing topography information.

        The topography data file can either have 
        two columns: profile position and topography,
        or three columns: X, Y, and topography, or Easting, Northing, topography
        
        In the topo text file, units of profile position (or northing and easting)
        and of the topography (or elevation) need to be in meters!

        Requires the velocity to be set.

        INPUT:
        topofile      file name for ASCII text topography information
        delimiter     how the entries are delimited (by comma, or by tab)
                      [default: ',']. To set tab: delimiter='\t'
        '''
        if self.velocity is None:
            print("First need to set velocity!")
            return
        # Store previous state for undo
        self.storePrevious()
        self.data_pretopo = self.data
        self.twtt_pretopo = self.twtt
        topoPos, topoVal, self.threeD = tools.prepTopo(topofile,delimiter,self.profilePos[0])
        self.data, self.twtt, self.maxTopo, self.minTopo = tools.correctTopo(self.data,
                                                                             velocity=self.velocity,
                                                                             profilePos=self.profilePos,
                                                                             topoPos=topoPos,
                                                                             topoVal=topoVal,
                                                                             twtt=self.twtt)
        # Put in history
        if delimiter is ',':
            histstr = "mygpr.topoCorrect('%s')" %(topofile)
        else:
            histstr = "mygpr.topoCorrect('%s',delimiter='\\t')" %(topofile)
        self.history.append(histstr)
        


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
        x,y,z = tools.prepVTK(self.profilePos,gpsmat,smooth,win_length,porder)        
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
        
        gridToVTK(outfile,XX,YY,ZZ, cellData ={'gpr': datarray})
 
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
        
    def storePrevious(self):
        '''
        Stores the current state of the profile and history in the 
        "previous" variable to be able to apply undo.
        '''
        self.previous["data"] = self.data
        self.previous["twtt"] = self.twtt
        self.previous["info"] = self.info
        self.previous["profilePos"] = self.profilePos
        self.previous["history"] = self.history
        self.previous["velocity"] = self.velocity
        self.previous["depth"] = self.depth
        self.previous["maxTopo"] = self.maxTopo
        self.previous["threeD"] = self.threeD
        self.previous["data_pretopo"] = self.data_pretopo
        self.previous["twtt_pretopo"] = self.twtt_pretopo
        



        
class gprpyCW(gprpyProfile):
    '''
    Ground penetrating radar data processing and visualization class
    for common midpoint or wide angle reflection and refraction data
 
    Inherits all of the gprpyProfile class functions but not all
    of these functions may be useful here. 
    '''
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
        super().__init__(None)
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


    def storePrevious(self):
        '''
        Stores the current state of the profile and history in the 
        "previous" variable to be able to apply undo.
        '''
        self.previous["data"] = self.data
        self.previous["twtt"] = self.twtt
        self.previous["info"] = self.info
        self.previous["profilePos"] = self.profilePos
        self.previous["history"] = self.history
        self.previous["dtype"] = self.dtype
        self.previous["hypStAmp"] = self.hypStAmp
        self.previous["linStAmp"] = self.linStAmp
        self.previous["lins"] = self.lins
        self.previous["hyps"] = self.hyps



            
    def importdata(self,filename,dtype):
        '''
        Loads .gpr (native GPRPy), .DT1 (Sensors and Software),
        .DZT (GSSI), .GPRhdr (ENVI standard BSQ), .rad (MALA)
        data files and populates all the gprpyProfile fields.

        INPUT: 
        filename  name of the .gpr, .DT1, .DZT, .GPRhdr, dat, rd3, 
                  or .rad file you want to import.
                  The header file name and the data file name 
                  have to be the same!
        dtype     data type. Either "CMP" or "WARR
        
        '''
        print(filename)
        super().importdata(filename)
        self.dtype = dtype
        self.vVals = None
        # Remove the history string from the super-class importing
        del self.history[-1]
        # Put what you did in history
        histstr = "mygpr.importdata('%s',dtype='%s')" %(filename,dtype)
        self.history.append(histstr)  


    # def setZeroTimeCW(self,newZeroTime):
    #     '''
    #     Deletes all data recorded before newZeroTime and 
    #     sets newZeroTime to zero.

    #     INPUT:
    #     newZeroTime     The new zero-time
    #     '''
    #     # Store previous state for undo
    #     self.storePrevious()
    #     # Find index of value that is nearest to newZeroTime
    #     zeroind = np.abs(self.twtt - newZeroTime).argmin()
    #     # Cut out everything before
    #     self.twtt = self.twtt[zeroind:] - newZeroTime
    #     #self.data = self.data[zeroind:,:]
    #     # Put what you did in history
    #     histstr = "mygpr.setZeroTime(%g)" %(newZeroTime)
    #     self.history.append(histstr)  
        
    
    def normalize(self):
        '''
        Divides each trace by its total energy to counteract the 
        loss of energy for wider antennae separations.
        '''
        # Store previous state for undo
        self.storePrevious()
        # Calculate norm of each trace and divide each trace by it
        self.data = np.divide(self.data,np.maximum(np.linalg.norm(self.data,axis=0),1e-8))
        print("normalized data set")
        # Put what you did in history
        histstr = "mygpr.normalize()"
        self.history.append(histstr) 


    def linStackedAmplitude(self,vmin=0.01,vmax=0.35,vint=0.01):
        '''
        Calculates the linear stacked amplitudes for each 
        travel time sample and the provided velocity range 
        by summing the pixels of the data that follow a line given 
        by the travel time zero offset and the velocity.

        INPUT:
        vmin       minimal velocity for which to calculate the 
                   stacked amplitude, in m/ns [default = 0.01 m/ns]
        vmax       maximum velocity for which to calculate the 
                   stacked amplitude, in m/ns [default = 0.35 m/ns]
        vint       velocity intervall, in m/ns [default = 0.01 m/ns]
        '''
        # Store previous state for undo
        self.storePrevious()
        self.vVals = np.arange(vmin,vmax+vint,vint)
        if self.dtype is "WARR":
            typefact = 1
        elif self.dtype is "CMP":
            typefact = 2
        self.linStAmp = tools.linStackedAmplitude(self.data,self.profilePos,self.twtt,self.vVals,self.twtt,typefact)
        print("calculated linear stacked amplitude")
        # Put what you did in history
        histstr = "mygpr.linStackedAmplitude(vmin=%g,vmax=%g,vint=%g)" %(vmin,vmax,vint)
        self.history.append(histstr)
        

    def hypStackedAmplitude(self,vmin=0.01,vmax=0.35,vint=0.01):
        '''
        Calculates the hyperbolic stacked amplitudes for each 
        travel time sample and the provided velocity range 
        by summing the pixels of the data that follow a hyperbola given 
        by the travel time apex and the velocity.

        INPUT:
        vmin       minimal velocity for which to calculate the 
                   stacked amplitude, in m/ns [default = 0.01 m/ns]
        vmax       maximum velocity for which to calculate the 
                   stacked amplitude, in m/ns [default = 0.35 m/ns]
        vint       velocity intervall, in m/ns [default = 0.01 m/ns]
        '''
        # Store previous state for undo
        self.storePrevious()
        self.vVals = np.arange(vmin,vmax+vint,vint)
        if self.dtype is "WARR":
            typefact = 1
        elif self.dtype is "CMP":
            typefact = 2
        self.hypStAmp = tools.hypStackedAmplitude(self.data,self.profilePos,self.twtt,self.vVals,self.twtt,typefact)
        print("calculated hyperbola stacked amplitude")
        # Put what you did in history
        histstr = "mygpr.hypStackedAmplitude(vmin=%g,vmax=%g,vint=%g)" %(vmin,vmax,vint)
        self.history.append(histstr)                  


    def addLin(self,zerotwtt,vel):
        '''
        Adds an observed line given by its zero-offset travel
        time and velocity to the list of lines.

        INPUT:
        zerotwtt     the zero-offset travel time (intercept)
                     of the observed line
        vel          the velocity (inverse slope) of the observed line
        '''
        # Store previous state for undo
        self.storePrevious()
        self.lins.append([zerotwtt,vel])
        # Put what you did in history
        histstr = "mygpr.addLin(zerotwtt=%g,vel=%g)" %(zerotwtt,vel)
        self.history.append(histstr)  
            

    def addHyp(self,zerotwtt,vel):
        '''
        Adds an observed hyperbola given by its apex 
        travel time and velocity to the list of lines.

        INPUT:
        zerotwtt     the apex travel time of the observed line
        vel          the velocity of the observed line
        '''
        # Store previous state for undo
        self.storePrevious()
        self.hyps.append([zerotwtt,vel])
        # Put what you did in history
        histstr = "mygpr.addHyp(zerotwtt=%g,vel=%g)" %(zerotwtt,vel)
        self.history.append(histstr)  

    def remLin(self):
        '''
        Removes the most recently added observed line from the list
        of observed lines
        '''
        # Store previous state for undo
        self.storePrevious()
        del self.lins[-1]
        # Put what you did in history
        histstr = "mygpr.remLin()" 
        self.history.append(histstr) 

    def remHyp(self):
        '''
        Removes the most recently added observed hyperbola from the list
        of observed lines
        '''
        # Store previous state for undo
        self.storePrevious()
        del self.hyps[-1]
        # Put what you did in history
        histstr = "mygpr.remHyp()" 
        self.history.append(histstr) 




    # This is a helper function
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
        
       

    def printCWFigure(self, figname, dpi=600, **kwargs):
        '''
        Creates a pdf of the common midpoint or wide angle reflection
        and refraction data.

        INPUT:
        figname      file name for the pdf
        dpi          dots per inch resolution [default: 600 dpi]
        color        "gray", or "bwr" for blue-white-red,
                     or any other Matplotlib color map [default: "gray"]
        contrast     Factor to increase contrast by reducing color range.
                     [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        xrng         x-axis range to show [default: None, meaning "everything"]
        showlnhp     show the observed lines and hyperbolae from the list
                     [default: False]
        '''

        
        contrast, color, yrng, xrng, showlnhp = self.prepCWFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        # Put what you did in history
        histstr = "mygpr.printCWFigure('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], showlnhp=%r, dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],showlnhp,dpi)
        self.history.append(histstr)


    def printStAmpFigure(self, figname, dpi=600, **kwargs):
        '''
        Creates a pdf of the the stacked amplitude results.

        INPUT:
        figname      file name for the pdf
        dpi          dots per inch resolution [default: 600 dpi]
        whichstamp   is this for the linear ("lin") or hyperbolic ("hyp") 
                     stacked amplitudes
        saturation   Factor to increase contrast by reducing color range.
                     [default = 1.0]
        yrng         y-axis range to show [default: None, meaning "everything"]
        vrng         velocities (x-axis) range to show 
                     [default: None, meaning "everything"]
        '''
        
        whichstamp, saturation, yrng, vrng = self.prepStAmpFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        histstr = "mygpr.printStAmpFigure('%s', whichstamp='%s', saturation=%g, yrng=[%g,%g], vrng=[%g,%g])" %(figname, whichstamp, saturation, yrng[0], yrng[1], vrng[0], vrng[1])
        self.history.append(histstr)
