from distutils.command import sdist
from turtle import fd
import tkinter.simpledialog as sd
import matplotlib as mpl
import numpy as np

class Visualization:

    def __init__(self):
        pass

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
            figrowsp = 1  # Define the variable figrowsp
            if self.picking:
                a.plot(self.picked[:,0],self.picked[:,1],'-x',color='yellow',linewidth=3*self.highfac) 
                a.plot(self.picked[:,0],self.picked[:,1],'-x',color='black',linewidth=2*self.highfac)                               

            # Allow for cursor coordinates being displayed        
            def moved(event):
                if event.xdata is not None and event.ydata is not None:
                    canvas.get_tk_widget().itemconfigure(tag, text="(x = %5.5g, y = %5.5g)" % (event.xdata, event.ydata))

            self.cursor_cid = canvas.mpl_connect('button_press_event', moved)
            tag = canvas.get_tk_widget().create_text(20, 20, text="", anchor="nw")

            canvas.get_tk_widget().grid(row=2,column=0,columnspan=figcolsp, rowspan=figrowsp, sticky='nsew')
            canvas.draw()
            

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

    def undo(self,proj):
        if self.picking:
            self.picked=self.picked[0:-1,:]
        else:
            proj.undo()   


    def setYrng(self, ylow, yhigh):
        self.prevyrng = self.yrng
        self.yrng = [ylow, yhigh]

    def setXrng(self, xlow, xhigh):
        self.xrng = [xlow, xhigh]

    def setVelRng(self):
        vmin = sd.askfloat("Input","Minimum velocity")
        if vmin is not None:
            self.vmin = vmin
        vmax = sd.askfloat("Input","Maximum velocity")
        if vmax is not None:
            self.vmax = vmax
        vint = sd.askfloat("Input","Velocity step size (interval)")
        if vint is not None:
            self.vint = vint    
            
    def setFullView(self,proj):
        self.xrng=[np.min(proj.profilePos),np.max(proj.profilePos)]
        self.yrng=[np.min(proj.twtt),np.max(proj.twtt)]