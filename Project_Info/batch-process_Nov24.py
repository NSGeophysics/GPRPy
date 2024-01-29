#############################################################
# Processing sequence proposed by Prof. J. Christian Dupuis 
# Geophysical Instrumentation Group Université Laval (GiGUL)
# -----------------------------------------------------------
# All of the raw data should be in the raw data folder
# The script does the rest 
#############################################################
from threading import _DummyThread
import gprpy.gprpy as gp
import matplotlib.pyplot as plt
import numpy as np
import os
import gprpy.makeDataCube as dc
import gprpy.interpSurface as ints
import datetime

################################################################
author = 'JCD'
analyst = 'WW' #Change Name of Analyst depending who did the GPR Analysis
Version =  0.1 
current_date_time = datetime.datetime.now()

# Initialise the data structure that holds the GPR data #######
mygpr = gp.gprpyProfile()
gpr_profile_for_interp = gp.gprpyProfile()
###############################################################

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

###############################################################
# Function that allows to generate the directories required
###############################################################

def setup_file_directories(dname): 
    try:
        os.makedirs(dname)
    except FileExistsError:
    # directory already exists
        pass
##############################################################

def reset_folder(dname):
    file_list = os.listdir(dname)
    for file in file_list:
        os.remove(dname+file)

# Try to setup the directory structure if it is required #####
setup_file_directories(proc_data_folder)
setup_file_directories(f_topo)
setup_file_directories(report_doc)
setup_file_directories(report_figures)
setup_file_directories(report_figures_interp)
##############################################################
# clean up folders if they had files in them
reset_folder(f_topo) 
reset_folder(proc_data_folder)
reset_folder(report_figures_interp)

########################################################
# Generate a list of profiles that we want to process ###
raw_grid = os.listdir(raw_data_folder) # Get the list of files in the raw folder
########################################################

def create_topo (xstart,xstop,ystart,ystop,zstart,zstop,raw_grid,f_topo):
# Generate the geometry information 
    num_step_y = len(raw_grid) # this is the number of profiles
    y,dy=np.linspace(ystart,ystop,num_step_y,endpoint=True, retstep=True)
    x = [xstart,xstop]
    z = [zstart,zstop]

    i = 0 
    for profile in raw_grid:
        profile_name,file_format = profile.split(".",1)
        topo = np.array([[x[0],y[i],z[0]],[x[1],y[i],z[1]]])
        np.savetxt(f_topo+'topo-'+profile_name+'.txt',topo,delimiter='\t',fmt='%.2f')
        print(topo)
        i=i+1
    return dy

##########################################################
# Create the topography for the arbitrary grid 
# used in acquisition
##########################################################
dy = create_topo (xstart,xstop,ystart,ystop,zstart,zstop,raw_grid,f_topo)


##########################################################
# Process the raw profiles that are in the raw_data_folder
# ######################################################## 
i = 0
for profile in raw_grid:
    mygpr.importdata(raw_data_folder+profile)
    mygpr.setZeroTime(zero_time)
    mygpr.adjProfile(0,30) # Adjust profile if distance wheel did not record properly
    profile_name,file_format = profile.split(".",1)
    plt.subplot(3,1,1)
    plot_title = ' ( x = %2.2f m )'%(i*xinc) + ' : ' + profile_name
    plt.title(plot_title)
    mygpr.showProfile(color='gray', contrast=4, yrng=[0,100], xrng=[0,xstop], asp=0.02) #Profile length is determined in Xrang from 0 to the max length of the profile (Line 45/ystop). Adjust the other ranges as needed.
    plt.ylabel('RAW')
    mygpr.remMeanTrace(ntraces)
    mygpr.dewow(dwowsamples)
    #mygpr.truncateY(20)
    mygpr.agcGain(agc)
    plt.subplot(3,1,2)
    mygpr.showProfile(color='gray', contrast=4, yrng=[0,100], xrng=[0,xstop], asp=0.02) 
    plt.ylabel('filtered')
    mygpr.setVelocity(v_migration)
    mygpr.fkMigration()
    plt.subplot(3,1,3)
    mygpr.showProfile(color='gray', contrast=4, yrng=[0,4], xrng=[0,xstop], asp=1)
    plt.ylabel('migrated')
    mygpr.topoCorrect(f_topo+'topo-'+profile_name+'.txt',delimiter='\t')
    mygpr.save(proc_data_folder + profile_name + 'processed' )
    #- Save the figure that document the processing steps ---------------------------
    plt.tight_layout() #No overlap with the axis labels
    plt.savefig(report_figures+profile_name+'.png',dpi=300)
    # -------------------------------------------------------------------------------
    i = i+1

########################################################
# Generate a list of profiles that we want to process ###
proc_grid = os.listdir(proc_data_folder) # Get the list of files in the raw folder
########################################################


###############################################################
# Prepare the images for interpretation
###############################################################
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

#####################################################################
# Write a report on the processing parameters and Metadata used to process the files
#####################################################################
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
    print("---------------------------------------------------------", file=f)
    print ("Removed mean traces using %i traces \n" %ntraces, file=f)
    print ("Data has agc gain curve applied with a window of %s \n" %agc, file=f)      ### Change these when using agc
    print ("Data is dewow'ed with a window size of %i \n"%dwowsamples, file=f)         ### Change these when using dwow
    print ("Data was migrated (f-k) using a constant velocity of %2.2f m/ns \n" % v_migration, file=f)
    print ("The processed data is in : %s \n" %proc_data_folder, file=f)
    print ("Figures for each profile as processed is in : %s \n" %report_figures, file=f)
    print ("3D data cube generated is called Cube.vts and is found in folder ./proc/", file=f)
    print("---------------------------------------------------------", file=f)


###################################################################
# Make the data cube with all of the profiles that were processed 
###################################################################
proc_grid = os.listdir(proc_data_folder)
datalist = list()
for profile in proc_grid:
    datalist.append(proc_data_folder+profile)

dc.makeDataCube(datalist,'Suckerville',nx=nx,ny=ny,nz=nz,smooth=(0.2,0.5,0.2),absvals=True) #Changenamefor fileoutput Normal for Gravses (0.2,0.5,0.2)
