import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gprIO_DT1
import gprpyTools as tools
import copy

class gprpy2d:
    def __init__(self,filename=None,desciption=None): #,profilerange=None):
        self.history = ["mygpr = gprpy.gprpy2d()"]

        # Initialize previous for undo
        self.previous = {
            "data" : None,
            "twtt" : None,
            "info" : None,
            "profilePos" : None,
            "history" : None}
        
        if filename is not None:
            self.importdata(filename)                 
        
    def importdata(self,filename):
        file_name, file_ext = os.path.splitext(filename)
        if file_ext==".DT1":
            self.data=gprIO_DT1.readdt1(filename)
            self.info=gprIO_DT1.readdt1Header(file_name + ".HD")
            
            self.profilePos = np.linspace(self.info["Start_pos"],
                                          self.info["Final_pos"],
                                          self.info["N_traces"])

            self.twtt = np.linspace(self.info["TZ_at_pt"],
                                    self.info["Total_time_window"],
                                    self.info["N_pts_per_trace"])
                
            # Put what you did in history
            histstr = "mygpr.importdata('%s')" %(filename)
            self.history.append(histstr)
            
        elif file_ext==".DZT":
            print("DZT Not yet implemented")
            
        elif file_ext==".gpr":
            #print("Not yet ready")
            ## Getting back the objects:
            with open(filename, 'rb') as f:
                data, info, profilePos, twtt, history = pickle.load(f)
            self.data = data
            self.info = info
            self.profilePos = profilePos
            self.twtt = twtt
            self.history = history
            
        else:
            print("Can only read dt1 or dzt files")

    def showHistory(self):
        for i in range(0,len(self.history)):
            print(self.history[i])

    def writeHistory(self,outfilename="myhistory.py"):
        with open(outfilename,"w") as outfile:
            outfile.write("import gprpy\n")
            for i in range(0,len(self.history)):
                outfile.write(self.history[i] + "\n")
                
    def undo(self):
        self.data = self.previous["data"]
        self.twtt = self.previous["twtt"]
        self.info = self.previous["info"]
        self.profilePos = self.previous["profilePos"]
        # Make sure to not keep deleting history
        # when applying undo several times. 
        histsav = copy.copy(self.previous["history"])
        del histsav[-1]
        self.history = histsav
        print("undo")
        

    def save(self,filename):
        # Saving the objects:
        # Want to force the file name .gpr
        file_name, file_ext = os.path.splitext(filename)
        if not(file_ext=='.gpr'):
            filename = filename + '.gpr'
        with open(filename, 'wb') as f:  
            pickle.dump([self.data, self.info, self.profilePos, self.twtt, self.history], f)
        print("Saved " + filename)
        # Add to history string
        histstr = "mygpr.save('%s')" %(filename)
        self.history.append(histstr)

        
    #def setRange(self, profilerange):
    #    # Only use this if the step size is not accurate
    #    self.profilerange=[min(profilerange),max(profilerange)]
    #    histstr = "mygpr.setRange([%f, %f])" %(min(profilerange),max(profilerange))
    #    self.history.append(histstr)


    
    # This is a helper function
    def prepTWTTfig(self, color="gray", contrast=1.0, timelim=None, profilelim=None):
        stdcont = np.argmax(abs(self.data)) 
        plt.imshow(self.data,cmap=color,extent=[min(self.profilePos),
                                                max(self.profilePos),
                                                max(self.twtt),
                                                min(self.twtt)],
                   aspect="auto",vmin=-stdcont/contrast, vmax=stdcont/contrast)
        if timelim is not None:
            plt.ylim(timelim)
            plt.gca().invert_yaxis()
        if profilelim is not None:
            plt.xlim(profilelim)
        #plt.gca().set_ylim([0,min(maxyval,max(proj.twtt))])
        #plt.gca().invert_yaxis()
        plt.gca().get_xaxis().set_visible(True)
        plt.gca().get_yaxis().set_visible(True)
        plt.gca().set_ylabel("two-way travel time [ns]")
        plt.gca().set_xlabel("profile position")
        plt.gca().xaxis.tick_top()
        plt.gca().xaxis.set_label_position('top')
        
        return contrast, color, timelim, profilelim
       
    
    def showTWTT(self, **kwargs):
        self.prepTWTTfig(**kwargs)
        plt.show(block=False)


    def printTWTT(self, figname, **kwargs):
        contrast, color, timelim, profilelim = self.prepTWTTfig(**kwargs)
        plt.savefig(figname, format='pdf')
        plt.close('all')
        # Put what you did in history
        histstr = "mygpr.printTWTT('%s', color='%s', contrast=%f, timelim=%s, profilelim=%s)" %(figname,color,contrast,timelim,profilelim)
        self.history.append(histstr)
        

    ####### Processing #######

    def timeZeroAdjust(self):
        # Store previous state for undo
        self.storePrevious()
        
        self.data = tools.timeZeroAdjust(self.data)
        
        # Put what you did in history
        histstr = "mygpr.timeZeroAdjust()"
        self.history.append(histstr)


    def dewow(self,window):
        # Store previous state for undo
        self.storePrevious()

        self.data = tools.dewow(self.data,window)

        # Put in history
        histstr = "mygpr.dewow(%d)" %(window)
        self.history.append(histstr)


    def remMeanTrace(self,ntraces):
        # Store previous state for undo
        self.storePrevious()

        self.data = tools.remMeanTrace(self.data,ntraces)
        
        # Put in history
        histstr = "mygpr.remMeanTrace(%d)" %(ntraces)
        self.history.append(histstr)


    def toDepth(self,velocity):
        # Store previous state for undo
        self.storePrevious()
        
        self.twtt = self.twtt * velocity/2.0

        # Put in history
        histstr = "mygpr.toDepth(%f)" %(velocity)
        self.history.append(histstr)
            


    def storePrevious(self):        
        self.previous["data"] = self.data
        self.previous["twtt"] = self.twtt
        self.previous["info"] = self.info
        self.previous["profilePos"] = self.profilePos
        self.previous["history"] = self.history
