import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import gprpy.gprpy as gp
import matplotlib.pyplot as plt
import numpy as np
import os
import gprpy.makeDataCube as dc
import gprpy.interpSurface as ints
import datetime

'''
# FOLDERS THAT CONTAIN DATA RAW AND PROCESSED DATA EVENTUALLY #
raw_data_folder = '../GPR_Data'  
proc_data_folder = 'proc/data_processed/'
f_topo = 'proc/data_topo/'
report_doc ='proc/report/'
report_figures = 'proc/report/figures/'
report_figures_interp = '../report/figures/interp/'
###############################################################

# parameters to setup for batch processing ####################
zero_time = 13.55 # time offset in ns
ntraces = 965 # number of traces to use in remove mean trace
v_migration = 0.07 # velocity of migration m/ns
dwowsamples = 1024  #Dewow Sample window size
agc = 1024   #agc gain window size
#  -----------------------------------------------------------#

# Information about the acquisition grid ######################
xstart = 0.0 # m 
xstop = 20.0 # m   # Add Length of the Profile / meters along profiling direction
ystart = 0.0 # m    
ystop = 7.0 # m   # This is the length of you x-axis or the axis that you move the GPR every x cm/m.
zstart = 0.0 # m
zstop = 0.0 # m
xinc = 0.5 # m Survey interval (The step distance you are moving the GPR) 
#-------------------------------------------------------------#

# Mesh parameters for the 3D cube #############################
nx = 965 # number of mesh in x   (Number of Traces)
ny = 15 # number of mesh in y   (Number of Profiles/Lines/Files)
nz = 1024 # number of mesh in z (Number of Samples)


https://stackoverflow.com/questions/51877124/how-to-select-a-directory-and-store-it-into-a-variable-in-tkinter
https://www.geeksforgeeks.org/python-tkinter-entry-widget/
https://stackoverflow.com/questions/17290717/using-file-path-from-gui-textbox-in-code-in-python
https://stackoverflow.com/questions/22769634/how-to-print-file-path-in-text-box-using-tkinter
https://www.geeksforgeeks.org/python-tkinter-toplevel-widget/
https://www.tutorialspoint.com/python/tk_toplevel.htm
'''
legal_fmt = ['DZT','dt1', 'DT1', 'hd', 'HD', 'dat', 'GPRhdr', 'rad', 'rd3', 'rd7','gpr']
class Info_Collect(tk.Toplevel):
    
    def __init__(self):
        super().__init__()
        self.title("Batch Process GPR Data")
        self.geometry('800x500')
        self.f_topo = '/proc/data_topo/'
        self.report_doc ='/proc/report/'
        self.report_figures = '/proc/report/figures/'
        self.report_figures_interp = '/proc/report/figures/interp/'

        raw_data_lb = ttk.Label(self, text = 'Raw data folder')
        self.raw_data = StringVar()
        raw_data_E = ttk.Entry(self, textvariable= self.raw_data)
        self.raw_data_btn = Button(self, text='Browse', command= self.setRaw)
        raw_data_lb.grid(row = 0, column= 0)
        raw_data_E.grid(row= 0, column= 1)
        self.raw_data_btn.grid(row=0, column= 2)

        proc_data_lb = ttk.Label(self, text = 'Processed Data destination')
        self.proc_data = StringVar()
        proc_data_E = ttk.Entry(self, textvariable= self.proc_data)
        self.proc_data_btn = Button(self, text='Browse', command= self.setProc)
        proc_data_lb.grid(row = 1, column= 0)
        proc_data_E.grid(row= 1, column= 1)
        self.proc_data_btn.grid(row= 1, column= 2)
        #----------------------------------------------
        ttk.Label(self, text = 'Batch Process Parameters ---------------').grid(row = 3, column=0)
        ttk.Label(self, text = 'Zero Time').grid(row=4, column=0)
        self.zt = StringVar()
        zt_E = ttk.Entry(self, textvariable= self.zt)
        zt_E.grid(row=4,column=1)

        ttk.Label(self, text = 'N Traces').grid(row=5, column=0)
        self.nt = StringVar()
        ntrace_E = ttk.Entry(self, textvariable= self.nt)
        ntrace_E.grid(row= 5, column=1)

        ttk.Label(self, text = 'Velocity of Migration').grid(row=6,column=0)
        self.vel = StringVar()
        vel_E = ttk.Entry(self, textvariable= self.vel)
        vel_E.grid(row=6, column=1)

        ttk.Label(self, text = 'Dewow Sample Window Size').grid(row=7,column=0)
        self.dwow = StringVar()
        dwow_E = ttk.Entry(self, textvariable= self.dwow)
        dwow_E.grid(row=7,column=1)

        ttk.Label(self, text = 'AGC Gain Window Size').grid(row=8,column=0)
        self.acg = StringVar()
        acg_E = ttk.Entry(self, textvariable= self.acg)
        acg_E.grid(row=8,column=1)
        #----------------------------------------------
        ttk.Label(self, text = 'Acquisition Grid Info ---------------').grid(row=10,column=0)
        ttk.Label(self, text = 'X - start(m)').grid(row=11,column=0)
        self.xstart = StringVar()
        xstart_E = ttk.Entry(self, textvariable= self.xstart)
        xstart_E.grid(row=11,column=1)

        ttk.Label(self, text = 'X - stop(m)').grid(row=11,column=2)
        self.xstop = StringVar()
        xstop_E = ttk.Entry(self, textvariable= self.xstop)
        xstop_E.grid(row=11,column=3)

        ttk.Label(self, text = 'Y - start(m)').grid(row=12,column=0)
        self.ystart = StringVar()
        ystart_E = ttk.Entry(self, textvariable= self.ystart)
        ystart_E.grid(row=12,column=1)

        ttk.Label(self, text = 'Y - stop(m)').grid(row=12,column=2)
        self.ystop = StringVar()
        ystop_E = ttk.Entry(self, textvariable= self.ystop)
        ystop_E.grid(row=12,column=3)

        ttk.Label(self, text = 'Z - start(m)').grid(row=13,column=0)
        self.zstart = StringVar()
        zstart_E = ttk.Entry(self, textvariable= self.zstart)
        zstart_E.grid(row=13,column=1)

        ttk.Label(self, text = 'Z - stop(m)').grid(row=13,column=2)
        self.zstop = StringVar()
        zstop_E = ttk.Entry(self, textvariable= self.zstop)
        zstop_E.grid(row=13,column=3)

        ttk.Label(self, text = 'Survey interval(x increment)').grid(row=14,column=0)
        self.xinc = StringVar()
        xinc_E = ttk.Entry(self, textvariable= self.xinc)
        xinc_E.grid(row=14,column=1)

        ttk.Label(self, text = 'Export format').grid(row=14,column=2)
        self.n = StringVar()
        expfmt = ttk.Combobox(self,width= 10, textvariable= self.n)
        expfmt['values'] = ('.png','.pdf')
        expfmt.grid(row=14, column=3)
        expfmt.current(0)

        #----------------------------------------------
        
        #dc_cb.invoke(self.enableDCparams)

        ttk.Label(self, text = 'Mesh Parameters for 3D Cube ').grid(row=16, column=0)
        ttk.Label(self, text = 'Number of Traces (nx)').grid(row=17, column=0)
        self.nx = StringVar()
        nx_E = ttk.Entry(self, textvariable= self.nx)
        nx_E.grid(row=17, column=1)
        

        ttk.Label(self, text = 'Number of Profiles (ny)').grid(row=17, column=2)
        self.ny = StringVar()
        ny_E = ttk.Entry(self, textvariable= self.ny)
        ny_E.grid(row=17, column=3)
        
        ttk.Label(self, text = 'Number of Samples (nz)').grid(row=18, column=0)
        self.nz = StringVar()
        nz_E = ttk.Entry(self,textvariable= self.nz)
        nz_E.grid(row=18, column=1)
        
        self.checked = tk.BooleanVar()
        self.checked.set(FALSE)
        dc_cb = ttk.Checkbutton(self, text= 'Create Data Cube', variable= self.checked)#, command=self.enableDCparams)
        dc_cb.grid(row=18,column=2)
        
        ttk.Button(self, text= 'RUN', command=self.runInterpolation).grid(row=18, column=3)
        #ttk.Text(self, height = 1, width = 3)
    def clickedBt(self):
        if(self.check_ip(self.nx)):  
            print('Empty in nx')
        else:
            print(self.nx.get())
        
    
    def check_ip(self, val):
        print('in the check')
        if val.get() == "":
            return 1

    def setRaw(self):
        fd = self.getDirectory()
        self.raw_data.set(fd)
        print(self.raw_data.get())
    
    def setProc(self):
        fd = self.getDirectory()
        self.proc_data.set(fd)
        print(self.proc_data.get())
    
    def getDirectory(self):
        fdest = filedialog.askdirectory()
        return fdest
        
    def setup_file_directories(self,dname): 
        #print(dname)
        try:
            os.makedirs(dname)
        except FileExistsError:
    # directory already exists
            pass
##############################################################

    def reset_folder(self,dname):
        file_list = os.listdir(dname)
        for file in file_list:
            os.remove(dname+file)

    def create_topo (self,raw_grid):
        ftopo_dest = self.proc_data.get() + self.f_topo
        print(ftopo_dest)
        # Generate the geometry information 
        num_step_y = len(raw_grid) # this is the number of profiles
        y,dy=np.linspace(float(self.ystart.get()),float(self.ystop.get()),num_step_y,endpoint=True, retstep=True)
        x = [float(self.xstart.get()),float(self.xstop.get())]
        z = [float(self.zstart.get()),float(self.zstop.get())]

        i = 0 
        for profile in raw_grid:
            profile_name,file_format= profile.split(".",1)
            topo = np.array([[x[0],y[i],z[0]],[x[1],y[i],z[1]]])
            np.savetxt(ftopo_dest+'topo-'+profile_name+'.txt',topo,delimiter='\t',fmt='%.2f')
            
            i=i+1
        return dy
    
    def procData(self, rawGrid,raw_data_folder, zero_time,xinc, xstop, ystop ,ntraces, dwowsamples, agc, v_migration,f_topo, proc_data_folder, report_figures, frmt):
        i = 0
        mygpr = gp.gprpyProfile()
        for profile in rawGrid:
            print(profile)
            profile_name,file_format= profile.split(".",1)
            if(file_format in legal_fmt):
                mygpr.importdata(raw_data_folder+'/'+profile)
                mygpr.setZeroTime(zero_time)
                mygpr.adjProfile(0,30) # Adjust profile if distance wheel did not record properly
                profile_name,file_format = profile.split(".",1)
                plt.subplot(3,1,1)
                plot_title = ' ( x = %2.2f m )'%(i*xinc) + ' : ' + profile_name
                plt.title(plot_title)
                mygpr.showProfile(color='gray', contrast=4, yrng=[0,18], xrng=[0,xstop], asp=0.5) #Profile length is determined in Xrang from 0 to the max length of the profile (Line 45/ystop). Adjust the other ranges as needed.
                plt.ylabel('RAW')
                mygpr.remMeanTrace(ntraces)
                mygpr.dewow(dwowsamples)
                #mygpr.truncateY(20)
                mygpr.agcGain(agc)
                plt.subplot(3,1,2)
                mygpr.showProfile(color='gray', contrast=4, yrng=[0,18], xrng=[0,xstop], asp=0.5) 
                plt.ylabel('filtered')
                mygpr.setVelocity(v_migration)
                mygpr.fkMigration()
                plt.subplot(3,1,3)
                mygpr.showProfile(color='gray', contrast=4, yrng=[0,0.7], xrng=[0,xstop], asp=10)
                plt.ylabel('migrated')
                mygpr.topoCorrect(f_topo+'topo-'+profile_name+'.txt',delimiter='\t')
                mygpr.save(proc_data_folder + profile_name + 'processed' )
                #- Save the figure that document the processing steps ---------------------------
                plt.tight_layout() #No overlap with the axis labels
                plt.savefig(report_figures+profile_name+'.png',dpi=300)
            # -------------------------------------------------------------------------------
            i = i+1

    ###############################################################
    # Prepare the images for interpretation
    ###############################################################
    def Interp_Prep(self, proc_grid,proc_data_folder, xinc,report_figures_interp, dy ):
        print("Prepping")
        gpr_profile_for_interp = gp.gprpyProfile()
        i = 0
        for profile in proc_grid:
            profile_name,file_format = profile.split(".",1)
            gpr_profile_for_interp.importdata(proc_data_folder+profile)
            #mygpr.showProfile(color='seismic', contrast=4, yrng=[0,50], xrng=[0,13], asp=0.05)
            gpr_profile_for_interp.showProfile(color='gray', contrast=4, yrng=[-2,0], xrng=[0,20], asp=1.2)
            plot_title = ' ( y = %2.2f m )'%(i*xinc) + ' : ' + profile_name
            plt.title(plot_title)
            plt.ylabel('Depth (m)')
            plt.savefig(report_figures_interp+profile_name+ '_' +'%2.2f-m'%(i*dy)+'.png',dpi=300)
            plt.close('all')
            i= i+1
    
    def makeCube(self,proc_data_folder,nx,ny,nz):
        print("making")
        proc_grid = os.listdir(proc_data_folder)
        datalist = list()
        for profile in proc_grid:
            datalist.append(proc_data_folder+profile)
        dc.makeDataCube(datalist,'Suckerville',nx=nx,ny=ny,nz=nz,smooth=(0.2,0.5,0.2),absvals=True) #Changenamefor fileoutput Normal for Gravses (0.2,0.5,0.2)



    def createReport(self):
        report_name = report_doc+'report-'+str(current_date_time.year) + '-' + str(current_date_time.month) + '-' +str(current_date_time.day) +'.txt' 
        with open(report_name, "a") as f:
            print("---------------------------------------------------------", file=f)
            print("Autogenerated Processing Report for IPIAFOA data", file=f)
            print("Using version %2.2f of the script developped by %s \n" %(Version,author), file=f)
            print("Processed on : %s \n" %current_date_time, file=f )
            print("Data Processed by %s \n" %(analyst), file=f)  #Change Name for the name of the GPR Analyst
            print("---------------------------------------------------------", file=f)
            print ("Found %i profiles in raw_data directory : % s \n" %(len(raw_grid),raw_data_folder), file=f)
            #print ("File format of the original files : %s \n" %file_format, file=f)
            print("---------------------------------------------------------", file=f)
            print ("GEOMETRY as defined in script", file=f)
            print("---------------------------------------------------------", file=f)
            print ("xmin : %2.2f m \n" %xstart, file=f)
            print ("xmax : %2.2f m \n" %xstop, file=f)
            print ("ymin : %2.2f m \n" %ystart, file=f)
            print ("ymax : %2.2f m \n" %ystop, file=f)
            print ("dy : %2.2f m \n" %dy, file=f)
            print ("zmax: %2.2f m \n" %zstart, file=f)
            print ("zmax: %2.2f m \n" %zstop, file=f)
            print ("---------------------------------------------------------", file=f)
            print ("Removed mean traces using %i traces \n" %ntraces, file=f)
            print ("Data has agc gain curve applied with a window of %s \n" %agc, file=f)      ### Change these when using agc
            print ("Data is dewow'ed with a window size of %i \n"%dwowsamples, file=f)         ### Change these when using dwow
            print ("Data was migrated (f-k) using a constant velocity of %2.2f m/ns \n" % v_migration, file=f)
            print ("The processed data is in : %s \n" %proc_data_folder, file=f)
            print ("Figures for each profile as processed is in : %s \n" %report_figures, file=f)
            print ("3D data cube generated is called Cube.vts and is found in folder ./proc/", file=f)
            print("---------------------------------------------------------", file=f)

  
    def runInterpolation(self):
        proc_dest = self.proc_data.get()
        proc_dest_dp = self.proc_data.get() + '/proc/data_processed/'
        self.setup_file_directories(proc_dest_dp )
        self.setup_file_directories(proc_dest + self.f_topo)
        self.setup_file_directories(proc_dest + self.report_doc)
        self.setup_file_directories(proc_dest + self.report_figures)
        self.setup_file_directories(proc_dest + self.report_figures_interp)
        rawGrid = os.listdir(self.raw_data.get())
        print(rawGrid)
        print("Before Create_topo")
        dy = self.create_topo(rawGrid)

        print("After proc data")
        self.procData(rawGrid,self.raw_data.get(),float(self.zt.get()),float(self.xinc.get()),float(self.xstop.get()),float(self.ystop.get()),float(self.nt.get()), float(self.dwow.get()),float(self.acg.get()), float(self.vel.get()), proc_dest+self.f_topo,proc_dest_dp, proc_dest+self.report_figures,self.n.get() )
        
        if(self.checked.get()):
            print("CUBING")
            proc_grid = os.listdir(proc_dest_dp) # Get the list of files in the data processed folder
            self.Interp_Prep(proc_grid, proc_dest_dp, float(self.xinc.get()), proc_dest+self.report_figures_interp, dy)
            self.makeCube(proc_dest_dp,int(self.nx.get()), int(self.ny.get()),int(self.nz.get()))
        
def create_extra():
    n_win = Info_Collect()
def main():
    
    root = tk.Tk()
    root.geometry('300x400')
    #n_win(root)
    tk.Button(root, text= 'click me', command= create_extra).grid(row = 0, column= 0)

    
    root.mainloop()

if __name__ == '__main__':
    main()