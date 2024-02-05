import sys
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

colsp=2
rightcol=9
halfwid=6
figrowsp=21+1
figcolsp=9


class DataManipulation:
    def __init__(self):
        pass

    def undo(self,proj):
        if self.picking:
            self.picked=self.picked[0:-1,:]
        else:
            proj.undo()                        
        
        
    def setYrng(self):
        ylow = sd.simpledialog.askfloat("Input","Min Y value",initialvalue=self.yrng[0])
        if ylow is not None:            
            yhigh = sd.simpledialog.askfloat("Input","Max Y value",initialvalue=self.yrng[1])
            if yhigh is not None:
                self.prevyrng=self.yrng
                self.yrng=[ylow,yhigh]
        

    def resetYrng(self,proj):
        # Only needed in undo, and only if what you want to
        # undo changed the y axis
        if ("setVelocity" in proj.history[-1]) or ("topoCorrect" in proj.history[-1]) and not self.picking: 
            self.yrng=self.prevyrng


    def setAspect(self):
        self.asp = sd.askfloat("Input","Plotting aspect ratio", initialvalue=self.asp)
        

    def setFullView(self,proj):    
        self.xrng=[np.min(proj.profilePos),np.max(proj.profilePos)]
        if proj.velocity is None:
            self.yrng=[np.min(proj.twtt),np.max(proj.twtt)]
        elif proj.maxTopo is None:
            self.yrng=[np.min(proj.depth),np.max(proj.depth)]
        else:
            self.yrng=[proj.minTopo-np.max(proj.depth),proj.maxTopo-np.min(proj.depth)]

            
    def toggleGrid(self):
        self.grid = not self.grid
            
            
    def setXrng(self):
        xlow = sd.askfloat("Input","Min X value",initialvalue=self.xrng[0])
        if xlow is not None:
            xhigh = sd.askfloat("Input","Max X value",initialvalue=self.xrng[1])
            if xhigh is not None:
                self.xrng=[xlow,xhigh]
        

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

                
    def setZeroTime(self,proj):
        newZeroTime = sd.askfloat("Input","New zero time")
        if newZeroTime is not None:
            proj.setZeroTime(newZeroTime=newZeroTime)
        
        
    def dewow(self,proj):
        window = sd.askinteger("Input","Dewow window width (number of samples)")
        if window is not None:
            proj.dewow(window=window)


    def smooth(self,proj):
        window = sd.askinteger("Input","Smoothing window width (number of samples)")
        if window is not None:
            proj.smooth(window=window)
            

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
        
    def cut(self,proj):
        minX = sd.simpledialog.askfloat("Input","Minimum profile position")
        if minX is not None:
            maxX = sd.simpledialog.askfloat("Input","Maximum profile position")
            if maxX is not None:
                proj.cut(minX,maxX)
            
    def setVelocity(self,proj):
        velocity =  sd.simpledialog.askfloat("Input","Radar wave velocity [m/ns]?")        
        if velocity is not None:
            proj.setVelocity(velocity)
            self.prevyrng=self.yrng
            self.yrng=[0,np.max(proj.depth)]

    def antennaSep(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Antenna Sep Error","You have to set the velocity first")
        proj.antennaSep()

            
    def fkMigration(self,proj):
        if proj.velocity is None:
            mesbox.showinfo("Migration Error","You have to set the velocity first")
        proj.fkMigration()


    def profileSmooth(self,proj):
        ntraces = sd.askinteger("Input","Smooth over how many traces?")
        if ntraces is not None:
            noversample = sd.askinteger("Input","Make how many copies of each trace?\nRecommended: Same as number of traces to be smoothed.")
            if noversample is not None:
                proj.profileSmooth(ntraces,noversample)
        
            
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




    def startPicking(self,proj,fig,a,canvas):
        self.picking = True
        self.picked = np.asmatrix(np.empty((0,2)))
        print("Picking mode on")
        def addPoint(event):
            self.picked = np.append(self.picked,np.asmatrix([event.xdata,event.ydata]),axis=0)
            self.plotProfileData(proj,fig=fig,a=a,canvas=canvas)
            print(self.picked)
        self.pick_cid = canvas.mpl_connect('button_press_event', addPoint)

            

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