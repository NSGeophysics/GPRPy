import gprpy.gprpy as gp
import numpy as np
from scipy.ndimage import zoom

def mergeProfiles(file1,file2,outfile,gapfill=0):
    '''
    Merges two GPR profiles by placing the second one at the end 
    of the first one. 

    Make sure you preprocessed them in GPRPy and save them to have the 
    correct starting and end times for the profile, or to both start at
    0 to just append the second profile at the end of the first profile.

    INPUT: 
    file1      File name (including path) of the first profile
    file2      File name (including path) of the second profile
    outfile    File name (including path) for the merged file
    gapfill    If there is a gap between the profiles, fill it with
               zeros (0) or NaN ('NaN')? [default: 0]
    '''

    # Load the two profiles
    profile1 = gp.gprpyProfile(file1)
    profile2 = gp.gprpyProfile(file2)

    # make sure starting and end times are the same
    assert (profile1.twtt[0]==profile2.twtt[0] and profile1.twtt[-1]==profile2.twtt[-1]), "\n\nUse GPRPy to cut the profiles to the same travel times\nCurrently: file 1 is %g ns to %g ns and file 2 is %g ns to %g ns \n" %(profile1.twtt[0],profile1.twtt[-1],profile2.twtt[0],profile2.twtt[-1])
    

    # If they don't have the same number of samples,
    # then we need to interpolate the data to make them fit
    if len(profile1.twtt) > len(profile2.twtt):
        zfac = len(profile1.twtt)/len(profile2.twtt)
        profile2.data = zoom(profile2.data,[zfac,1])
            
    elif len(profile1.twtt) < len(profile2.twtt):
        zfac = len(profile2.twtt)/len(profile1.twtt)
        profile1.data = zoom(profile1.data,[zfac,1])
        profile1.twtt = profile2.twtt

    # If they don't have the same along-profile sampling,
    # need to interpolate the data such that it makes sense:
    if np.diff(profile1.profilePos)[3] < np.diff(profile2.profilePos)[3]:      
        zfac = np.diff(profile2.profilePos)[3]/np.diff(profile1.profilePos)[3]
        profile2.data = zoom(profile2.data,[1,zfac])
        profile2.profilePos=zoom(profile2.profilePos,zfac)       

    elif np.diff(profile1.profilePos)[3] > np.diff(profile2.profilePos)[3]:
        zfac = np.diff(profile1.profilePos)[3]/np.diff(profile2.profilePos)[3]
        profile1.data = zoom(profile1.data,[1,zfac])
        profile1.profilePos=zoom(profile1.profilePos,zfac)
        
    # Now concatenate the profile positions
    # In case someone didn't adjust their profile but just tries to merge them:
    if abs(profile2.profilePos[0]) < 1e-5:
        profile2.profilePos = profile2.profilePos + profile1.profilePos[-1]+np.diff(profile2.profilePos)[1]    
        
    # Otherwise they probably know what they are doing
    # If there is a gap, create an array with zeros or NaNs
    dx=np.diff(profile2.profilePos)[0]
    if profile2.profilePos[0] - profile1.profilePos[-1] < dx:
        nfill = int(np.round((profile2.profilePos[0] -
                              profile1.profilePos[-1])/dx))
        posfill = np.arange(0,nfill)*dx + profile1.profilePos[-1] + dx
        datfill = np.empty(((profile2.data).shape[0],nfill))
        if gapfill == 0:
            datfill.fill(0)
        else:
            datfill.fill(np.NaN)
        #datfill = np.zeros(((profile2.data).shape[0],nfill))
        profile2.profilePos=np.append(posfill,profile2.profilePos)
        profile2.data = np.hstack((datfill,profile2.data))

    # Append profile positions    
    profile1.profilePos = np.append(profile1.profilePos,profile2.profilePos)
    
    # Now merge them into profile 1      
    profile1.data = np.asmatrix(np.hstack((profile1.data,profile2.data)))
    
    # Set history to shortest possible:
    profile1.history = ["mygpr = gp.gprpyProfile()", "mygpr.importdata('%s.gpr')" %(outfile)]

    profile1.info="Merged"

    # Save the result in a .gpr file
    profile1.save(outfile)
