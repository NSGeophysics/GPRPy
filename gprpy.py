import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gprIO_DT1

class gprpy2d:
    def __init__(self,filename=None,desciption=None): #,profilerange=None):
        self.history = ["mygpr = gprpy.gprpy2d()"]
        
        if filename is not None:
            self.importdata(filename)

        #if profilerange is not None:
        #   self.setRange(profilerange)
        
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
            histstr = "mygpr.importdata(" + filename + ")"
            self.history.append(histstr)
            
        elif file_ext==".DZT":
            pass
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


    #def save(filename)
    ## Saving the objects:
    #with open('objs.pkl', 'w') as f:  # Python 3: open(..., 'wb')
    #    pickle.dump([obj0, obj1, obj2], f)
    ## Getting back the objects:
    #with open('objs.pkl') as f:  # Python 3: open(..., 'rb')
    #    obj0, obj1, obj2 = pickle.load(f)

    #def setRange(self, profilerange):
    #    # Only use this if the step size is not accurate
    #    self.profilerange=[min(profilerange),max(profilerange)]
    #    histstr = "mygpr.setRange([%f, %f])" %(min(profilerange),max(profilerange))
    #    self.history.append(histstr)
            
    def showTWTT(self, color="gray", timelim=None, profilelim=None):
        plt.imshow(self.data,cmap=color,extent=[min(self.profilePos),
                                                max(self.profilePos),
                                                max(self.twtt),
                                                min(self.twtt)], aspect="auto")

        if timelim is not None:
            plt.ylim(timelim)
            plt.gca().invert_yaxis()

        if profilelim is not None:
            plt.xlim(profilelim)

        plt.show()
        
