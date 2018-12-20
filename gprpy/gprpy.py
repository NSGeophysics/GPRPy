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
    print("No fk migration in public version")
import copy
import scipy.interpolate as interp
from pyevtk.hl import gridToVTK

class gprpy2d:
    def __init__(self,filename=None):
        self.history = ["mygpr = gp.gprpy2d()"]

        # Initialize previous for undo
        self.previous = {}
        
        if filename is not None:
            self.importdata(filename)                 
        
    def importdata(self,filename):
        '''

        Loads .gpr (native GPRPy), .DT1 (Sensors and Software),
        .DZT (GSSI), .GPRhdr (ENVI standard BSQ), .rad (MALA)
        data files and populates all the gprpy2d fields.

        INPUT: 

        filename  name of the .gpr, .DT1, .DZT, .GPRhdr, or .rad file
                  you want to import.
                  The header file name and the data file name 
                  have to be the same!

        Last modified by plattner-at-alumni.ethz.ch, 12/18/2018

        '''
        
        file_name, file_ext = os.path.splitext(filename)
        
        if file_ext==".DT1" or file_ext==".HD":
            self.data=gprIO_DT1.readdt1(filename)
            self.info=gprIO_DT1.readdt1Header(file_name + ".HD")
            
            self.profilePos = np.linspace(self.info["Start_pos"],
                                          self.info["Final_pos"],
                                          self.info["N_traces"])

            self.twtt = np.linspace(self.info["TZ_at_pt"],
                                    self.info["Total_time_window"],
                                    self.info["N_pts_per_trace"])

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

            self.profilePos = self.info["startposition"]+np.linspace(0.0,
                                                                     self.data.shape[1]/self.info["scpmeter"],
                                                                     self.data.shape[1])
            
            self.twtt = np.linspace(0,self.info["nanosecptrace"],self.info["sptrace"])

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
                data, info, profilePos, twtt, history, velocity, depth, maxTopo, minTopo, threeD, data_pretopo, twtt_pretopo = pickle.load(f)
            self.data = data
            self.info = info
            self.profilePos = profilePos
            self.twtt = twtt
            self.history = history
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
            print("Can only read dt1 or dzt files")

    def showHistory(self):
        for i in range(0,len(self.history)):
            print(self.history[i])

    def writeHistory(self,outfilename="myhistory.py"):
        with open(outfilename,"w") as outfile:
            outfile.write("# Automatically generated by GPRPy\nimport gprpy.gprpy as gp\n")
            for i in range(0,len(self.history)):
                outfile.write(self.history[i] + "\n")
                
    def undo(self):
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
        # Saving the objects:
        # Want to force the file name .gpr
        file_name, file_ext = os.path.splitext(filename)
        if not(file_ext=='.gpr'):
            filename = filename + '.gpr'
        with open(filename, 'wb') as f:  
            pickle.dump([self.data, self.info, self.profilePos, self.twtt,
                         self.history, self.velocity, self.depth, self.maxTopo,
                         self.minTopo, self.threeD, self.data_pretopo,
                         self.twtt_pretopo], f)
        print("Saved " + filename)
        # Add to history string
        histstr = "mygpr.save('%s')" %(filename)
        self.history.append(histstr)

    
    # This is a helper function
    def prepProfileFig(self, color="gray", contrast=1.0, yrng=None, xrng=None, asp=None):
        dx=self.profilePos[1]-self.profilePos[0]
        dt=self.twtt[1]-self.twtt[0]
        stdcont = np.nanmax(np.abs(self.data)[:])       
        
        if self.velocity is None:
            plt.imshow(self.data,cmap=color,extent=[min(self.profilePos)-dx/2.0,
                                                    max(self.profilePos)+dx/2.0,
                                                    max(self.twtt)+dt/2.0,
                                                    min(self.twtt)-dt/2.0],
                       aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)
            plt.gca().set_ylabel("two-way travel time [ns]")
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
        
        return contrast, color, yrng, xrng, asp
       
    
    def showProfile(self, **kwargs):
        self.prepProfileFig(**kwargs)
        plt.show(block=False)


    def printProfile(self, figname, dpi=600, **kwargs):
        contrast, color, yrng, xrng, asp = self.prepProfileFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        # Put what you did in history
        if asp is None:
            histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],dpi)
        else:
            histstr = "mygpr.printProfile('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], asp=%g, dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],asp,dpi)
        self.history.append(histstr)
        

    ####### Processing #######

    def setRange(self,minPos,maxPos):
        # Adjust the length of the profile, in case the trigger wheel is not
        # Calibrated
        # Store previous state for undo
        self.storePrevious()
        self.profilePos=np.linspace(minPos,maxPos,np.size(self.profilePos))
        histstr = "mygpr.setRange(%g,%g)" %(minPos,maxPos)
        self.history.append(histstr)


    def flipProfile(self):
        # Flips the profile left to right (start to end)
        self.storePrevious()
        self.data=np.flip(self.data,1)
        histstr = "mygpr.flipProfile()"
        self.history.append(histstr)
        

    def alignTraces(self):
        # Store previous state for undo
        self.storePrevious()        
        self.data = tools.alignTraces(self.data)      
        # Put what you did in history
        histstr = "mygpr.alignTraces()"
        self.history.append(histstr)


    def adjProfile(self,minPos,maxPos):
        # Store previous state for undo
        self.storePrevious()
        # set new profile positions
        self.profilePos = np.linspace(minPos,maxPos,len(self.profilePos))       
        # Put what you did in history
        histstr = "mygpr.adjProfile(%g,%g)" %(minPos,maxPos)
        self.history.append(histstr)   

        
    def setZeroTime(self,newZeroTime):
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
        # Store previous state for undo
        self.storePrevious()
        self.data = tools.dewow(self.data,window)
        # Put in history
        histstr = "mygpr.dewow(%d)" %(window)
        self.history.append(histstr)


    def smooth(self,window):
        # Store previous state for undo
        self.storePrevious()
        self.data = tools.smooth(self.data,window)
        # Put in history
        histstr = "mygpr.smooth(%d)" %(window)
        self.history.append(histstr)

        
    def remMeanTrace(self,ntraces):
        # Store previous state for undo
        self.storePrevious()
        # apply
        self.data = tools.remMeanTrace(self.data,ntraces)        
        # Put in history
        histstr = "mygpr.remMeanTrace(%d)" %(ntraces)
        self.history.append(histstr)


    def profileSmooth(self,ntraces,noversample):
        # Store previous state for undo
        self.storePrevious()
        self.data,self.profilePos = tools.profileSmooth(self.data,self.profilePos,
                                                        ntraces,noversample)
        # Put in history
        histstr = "mygpr.profileSmooth(%d,%d)" %(ntraces,noversample)
        self.history.append(histstr)

        
    def tpowGain(self,power=0.0):
        # Store previous state for undo
        self.storePrevious()
        # apply tpowGain
        self.data = tools.tpowGain(self.data,self.twtt,power)
        # Put in history
        histstr = "mygpr.tpowGain(%g)" %(power)
        self.history.append(histstr)

    def agcGain(self,window=10):
        # Store previous state for undo
        self.storePrevious()
        # apply agcGain
        self.data = tools.agcGain(self.data,window)
        # Put in history
        histstr = "mygpr.agcGain(%d)" %(float(window))
        self.history.append(histstr)
        

    def setVelocity(self,velocity):
        # Store previous state for undo
        self.storePrevious()

        self.velocity = velocity
        self.depth = self.twtt * velocity/2.0

        # Put in history
        histstr = "mygpr.setVelocity(%g)" %(velocity)
        self.history.append(histstr)


    def fkMigration(self, **kwargs):
        # Store previous state for undo
        self.storePrevious()
        # apply migration
        dt=self.twtt[1]-self.twtt[0]
        dx=self.profilePos[1]-self.profilePos[0]
        # fkmig sets x profile to start at zero but resamples
        self.data,self.twtt,migProfilePos=mig_fk.fkmig(self.data,dt,dx,self.velocity)
        self.profilePos = migProfilePos + self.profilePos[0]
        
        # Put in history
        histstr = "mygpr.fkMigration()"
        self.history.append(histstr)
        
        
    def truncateY(self,maxY):
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
        


    def exportVTK(self,outfile,gpsinfo,thickness=0,delimiter=',',aspect=1.0,smooth=True, win_length=51, porder=3):
        # If gpsmat is a filename, we first need to load the file:
        if type(gpsinfo) is str:
            gpsmat = np.loadtxt(gpsinfo,delimiter=delimiter)
        else:
            gpsmat = gpsinfo
            
        # First get the x,y,z positions of our data points
        x,y,z = tools.prepVTK(self.profilePos,gpsmat,delimiter,smooth,win_length,porder)        
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
        



        
class gprpyCW(gprpy2d):
    def __init__(self,filename=None,dtype=None):
        '''
        dtype needs to be either "CMP" or "WARR"
        '''
        # Inheriting the initializer from the gprpy2d class
        super().__init__(None)
        self.history = ["mygpr = gp.gprpyCW()"]
        # Initialize previous for undo
        self.previous = {}
        self.dtype = dtype
        # Semblance plots
        self.linSemb = None
        self.hypSemb = None
        # Picked lines and hyperbolae
        self.lins = list()
        self.hyps = list()
        
        if (filename is not None) and (dtype is not None):
            self.importdata(filename,dtype)


    def storePrevious(self):        
        self.previous["data"] = self.data
        self.previous["twtt"] = self.twtt
        self.previous["info"] = self.info
        self.previous["profilePos"] = self.profilePos
        self.previous["history"] = self.history
        self.previous["dtype"] = self.dtype
        self.previous["hypSemb"] = self.hypSemb
        self.previous["linSemb"] = self.linSemb
        self.previous["lins"] = self.lins
        self.previous["hyps"] = self.hyps



            
    def importdata(self,filename,dtype):
        print(filename)
        super().importdata(filename)
        self.dtype = dtype
        self.vVals = None
        # Remove the history string from the super-class importing
        del self.history[-1]
        # Put what you did in history
        histstr = "mygpr.importdata('%s',dtype='%s')" %(filename,dtype)
        self.history.append(histstr)  


    def setZeroTimeCW(self,newZeroTime):
        # Store previous state for undo
        self.storePrevious()
        # Find index of value that is nearest to newZeroTime
        zeroind = np.abs(self.twtt - newZeroTime).argmin()
        # Cut out everything before
        self.twtt = self.twtt[zeroind:] - newZeroTime
        #self.data = self.data[zeroind:,:]
        # Put what you did in history
        histstr = "mygpr.setZeroTime(%g)" %(newZeroTime)
        self.history.append(histstr)  
        
    
    def normalize(self):
        # Store previous state for undo
        self.storePrevious()
        # Calculate norm of each trace and divide each trace by it
        self.data = np.divide(self.data,np.maximum(np.linalg.norm(self.data,axis=0),1e-8))
        print("normalized data set")
        # Put what you did in history
        histstr = "mygpr.normalize()"
        self.history.append(histstr) 


    def linSemblance(self,vmin=0.01,vmax=0.35,vint=0.01):
        # Store previous state for undo
        self.storePrevious()
        self.vVals = np.arange(vmin,vmax+vint,vint)
        if self.dtype is "WARR":
            typefact = 1
        elif self.dtype is "CMP":
            typefact = 2
        self.linSemb = tools.linSemblance(self.data,self.profilePos,self.twtt,self.vVals,self.twtt,typefact)
        print("calculated linear semblance")
        # Put what you did in history
        histstr = "mygpr.linSemblance(vmin=%g,vmax=%g,vint=%g)" %(vmin,vmax,vint)
        self.history.append(histstr)
        

    def hypSemblance(self,vmin=0.01,vmax=0.35,vint=0.01):
        # Store previous state for undo
        self.storePrevious()
        self.vVals = np.arange(vmin,vmax+vint,vint)
        if self.dtype is "WARR":
            typefact = 1
        elif self.dtype is "CMP":
            typefact = 2
        self.hypSemb = tools.hypSemblance(self.data,self.profilePos,self.twtt,self.vVals,self.twtt,typefact)
        print("calculated hyperbola semblance")
        # Put what you did in history
        histstr = "mygpr.hypSemblance(vmin=%g,vmax=%g,vint=%g)" %(vmin,vmax,vint)
        self.history.append(histstr)                  


    def addLin(self,zerotwtt,vel):
        # Store previous state for undo
        self.storePrevious()
        self.lins.append([zerotwtt,vel])
        # Put what you did in history
        histstr = "mygpr.addLin(zerotwtt=%g,vel=%g)" %(zerotwtt,vel)
        self.history.append(histstr)  
            

    def addHyp(self,zerotwtt,vel):
        # Store previous state for undo
        self.storePrevious()
        self.hyps.append([zerotwtt,vel])
        # Put what you did in history
        histstr = "mygpr.addHyp(zerotwtt=%g,vel=%g)" %(zerotwtt,vel)
        self.history.append(histstr)  

    def remLin(self):
        # Store previous state for undo
        self.storePrevious()
        del self.lins[-1]
        # Put what you did in history
        histstr = "mygpr.remLin()" 
        self.history.append(histstr) 

    def remHyp(self):
        # Store previous state for undo
        self.storePrevious()
        del self.hyps[-1]
        # Put what you did in history
        histstr = "mygpr.remHyp()" 
        self.history.append(histstr) 




    # This is a helper function
    def prepCWFig(self, contrast=1.0, color="gray", yrng=None, xrng=None, showlnhp=False):
        dx=self.profilePos[1]-self.profilePos[0]
        dt=self.twtt[1]-self.twtt[0]
        stdcont = np.nanmax(np.abs(self.data)[:])       
        
        plt.imshow(self.data,cmap=color,extent=[min(self.profilePos)-dx/2.0,
                                                max(self.profilePos)+dx/2.0,
                                                max(self.twtt)+dt/2.0,
                                                min(self.twtt)-dt/2.0],
                   aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)
        plt.gca().set_ylabel("two-way travel time [ns]")
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



    def prepSembFig(self, whichsemb="lin", saturation=1.0, yrng=None, vrng=None, sembrep="lin"):
        dt=self.twtt[1]-self.twtt[0]
        dv=self.vVals[1]-self.vVals[0]
        if whichsemb == "lin":
            semb = self.linSemb
            title = "linear semblance"
        elif whichsemb == "hyp":
            semb = self.hypSemb
            title = "hyperbolic semblance"
        else:
            semb = None
            
        if semb is not None:
            if sembrep == "lin":
                stdcont = np.nanmax(np.abs(semb)[:])
                plt.imshow(np.flipud(np.abs(semb)), cmap='inferno',
                           extent=[np.min(self.vVals)-dv/2.0, np.max(self.vVals)+dv/2.0,
                                   np.min(self.twtt)-dt/2.0,  np.max(self.twtt)+dt/2.0],
                           aspect='auto',
                           vmin=0, vmax=stdcont/saturation)
            elif self.sembrep.get() == "log":
                stdcont = np.nanmax(np.log(np.abs(semb))[:])
                plt.imshow(np.flipud(np.log(np.abs(semb))), cmap='inferno',
                           extent=[np.min(self.vVals)-dv/2.0, np.max(self.vVals)+dv/2.0,
                                   np.min(self.twtt)-dt/2.0,  np.max(self.twtt)+dt/2.0],
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
            plt.gca().set_ylabel("two-way travel time [ns]")
            #plt.gca().invert_yaxis()
            plt.gca().set_title(title)
            plt.gca().get_xaxis().set_visible(True)
            plt.gca().get_yaxis().set_visible(True)
            plt.gca().get_xaxis().set_ticks_position('both')
            plt.gca().get_yaxis().set_ticks_position('both')
                                
        return whichsemb, saturation, yrng, vrng, sembrep

    
    
    def showCWFig(self, **kwargs):
        self.prepCWFig(**kwargs)
        plt.show(block=False)


    def showSembFig(self, **kwargs):        
        self.prepSembFig(**kwargs)
        plt.show(block=False)
        
       

    def printCWFigure(self, figname, dpi=600, **kwargs):
        contrast, color, yrng, xrng, showlnhp = self.prepCWFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        # Put what you did in history
        histstr = "mygpr.printCWFigure('%s', color='%s', contrast=%g, yrng=[%g,%g], xrng=[%g,%g], showlnhp=%r, dpi=%d)" %(figname,color,contrast,yrng[0],yrng[1],xrng[0],xrng[1],showlnhp,dpi)
        self.history.append(histstr)


    def printSembFigure(self, figname, dpi=600, **kwargs):
        whichsemb, saturation, yrng, vrng, sembrep = self.prepSembFig(**kwargs)
        plt.savefig(figname, format='pdf', dpi=dpi)
        plt.close('all')
        histstr = "mygpr.printSembFigure('%s', whichsemb='%s', saturation=%g, yrng=[%g,%g], vrng=[%g,%g], sembrep='%s')" %(figname, whichsemb, saturation, yrng[0], yrng[1], vrng[0], vrng[1], sembrep)
        self.history.append(histstr)
