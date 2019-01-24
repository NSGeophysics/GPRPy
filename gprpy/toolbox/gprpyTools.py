import numpy as np
import scipy as sp
import numpy.matlib as matlib
import scipy.interpolate as interp
import scipy.signal as signal
# For progress bar
import time
from tqdm import tqdm


def alignTraces(data):
    '''
    Aligns the traces in the profile such that their maximum 
    amplitudes align at the average two-way travel time of the 
    maximum amplitudes 

    INPUT:
    data       data matrix whose columns contain the traces

    OUTPUT:
    newdata    data matrix with aligned traces
    '''
    
    maxlen = data.shape[0]
    newdata = np.asmatrix(np.zeros(data.shape))    
    # Go through all traces to find maximum spike
    maxind = np.zeros(data.shape[1], dtype=int)
    for tr in range(0,data.shape[1]):
        maxind[tr] = int(np.argmax(np.abs(data[:,tr])))
    # Find the mean spike point
    meanind = int(np.round(np.mean(maxind)))
    # Shift all traces. If max index is smaller than
    # mean index, then prepend zeros, otherwise append
    for tr in range(0,data.shape[1]):
        if meanind > maxind[tr]:
            differ = int(meanind - maxind[tr])
            newdata[:,tr] = np.vstack([np.zeros((differ,1)), data[0:(maxlen-differ),tr]])
        elif meanind < maxind[tr]:
            differ = maxind[tr] - meanind
            newdata[:,tr] = np.vstack([data[differ:maxlen,tr], np.zeros((differ,1))])
        else:
            newdata[:,tr] = data[:,tr]
    return newdata



def dewow(data,window):
    '''
    Subtracts from each sample along each trace an 
    along-time moving average.

    Can be used as a low-cut filter.

    INPUT:
    data       data matrix whose columns contain the traces 
    window     length of moving average window 
               [in "number of samples"]

    OUTPUT:
    newdata    data matrix after dewow
    '''
    totsamps = data.shape[0]
    # If the window is larger or equal to the number of samples,
    # then we can do a much faster dewow
    if (window >= totsamps):
        newdata = data-np.matrix.mean(data,0)            
    else:
        newdata = np.asmatrix(np.zeros(data.shape))
        halfwid = int(np.ceil(window/2.0))
        
        # For the first few samples, it will always be the same
        avgsmp=np.matrix.mean(data[0:halfwid+1,:],0)
        newdata[0:halfwid+1,:] = data[0:halfwid+1,:]-avgsmp

        # for each sample in the middle
        for smp in tqdm(range(halfwid,totsamps-halfwid+1)):
            winstart = int(smp - halfwid)
            winend = int(smp + halfwid)
            avgsmp = np.matrix.mean(data[winstart:winend+1,:],0)
            newdata[smp,:] = data[smp,:]-avgsmp

        # For the last few samples, it will always be the same
        avgsmp = np.matrix.mean(data[totsamps-halfwid:totsamps+1,:],0)
        newdata[totsamps-halfwid:totsamps+1,:] = data[totsamps-halfwid:totsamps+1,:]-avgsmp
        
    print('done with dewow')
    return newdata



def smooth(data,window):
    '''
    Replaces each sample along each trace with an 
    along-time moving average.

    Can be used as high-cut filter.

    INPUT:
    data      data matrix whose columns contain the traces 
    window    length of moving average window
              [in "number of samples"]

    OUTPUT:
    newdata   data matrix after applying smoothing
    '''
    totsamps = data.shape[0]
    # If the window is larger or equal to the number of samples,
    # then we can do a much faster dewow
    if (window >= totsamps):
        newdata = np.matrix.mean(data,0)
    elif window == 1:
        newdata = data
    elif window == 0:
        newdata = data
    else:
        newdata = np.asmatrix(np.zeros(data.shape))
        halfwid = int(np.ceil(window/2.0))
        
        # For the first few samples, it will always be the same
        newdata[0:halfwid+1,:] = np.matrix.mean(data[0:halfwid+1,:],0)

        # for each sample in the middle
        for smp in tqdm(range(halfwid,totsamps-halfwid+1)):
            winstart = int(smp - halfwid)
            winend = int(smp + halfwid)
            newdata[smp,:] = np.matrix.mean(data[winstart:winend+1,:],0)

        # For the last few samples, it will always be the same
        newdata[totsamps-halfwid:totsamps+1,:] = np.matrix.mean(data[totsamps-halfwid:totsamps+1,:],0)
        
    print('done with smoothing')
    return newdata



def remMeanTrace(data,ntraces):
    '''
    Subtracts from each trace the average trace over
    a moving average window.

    Can be used to remove horizontal arrivals, 
    such as the airwave.

    INPUT:
    data       data matrix whose columns contain the traces 
    ntraces    window width; over how many traces 
               to take the moving average.

    OUTPUT:
    newdata    data matrix after subtracting average traces
    '''

    data=np.asmatrix(data)
    tottraces = data.shape[1]
    # For ridiculous ntraces values, just remove the entire average
    if ntraces >= tottraces:
        newdata=data-np.matrix.mean(data,1) 
    else: 
        newdata = np.asmatrix(np.zeros(data.shape))    
        halfwid = int(np.ceil(ntraces/2.0))
        
        # First few traces, that all have the same average
        avgtr=np.matrix.mean(data[:,0:halfwid+1],1)
        newdata[:,0:halfwid+1] = data[:,0:halfwid+1]-avgtr
        
        # For each trace in the middle
        for tr in tqdm(range(halfwid,tottraces-halfwid+1)):   
            winstart = int(tr - halfwid)
            winend = int(tr + halfwid)
            avgtr=np.matrix.mean(data[:,winstart:winend+1],1)                
            newdata[:,tr] = data[:,tr] - avgtr

        # Last few traces again have the same average    
        avgtr=np.matrix.mean(data[:,tottraces-halfwid:tottraces+1],1)
        newdata[:,tottraces-halfwid:tottraces+1] = data[:,tottraces-halfwid:tottraces+1]-avgtr

    print('done with removing mean trace')
    return newdata



def profileSmooth(data,profilePos,ntraces=1,noversample=1):
    '''
    First creates copies of each trace and appends the copies 
    next to each trace, then replaces each trace with the 
    average trace over a moving average window.

    Can be used to smooth-out noisy reflectors appearing 
    in neighboring traces, or simply to increase the along-profile 
    resolution by interpolating between the traces.

    INPUT:
    data            data matrix whose columns contain the traces 
    profilePos      profile coordinates for the traces in data
    ntraces         window width [in "number of samples"]; 
                    over how many traces to take the moving average. 
    noversample     how many copies of each trace

    OUTPUT:
    newdata         data matrix after along-profile smoothing 
    newProfilePos   profile coordinates for output data matrix
    '''
    # New profile positions
    newProfilePos = np.linspace(profilePos[0],
                                profilePos[-1],
                                noversample*len(profilePos))
    # First oversample the data
    data = np.asmatrix(np.repeat(data,noversample,1))
    tottraces = data.shape[1]
    if ntraces == 1:
        newdata = data
    elif ntraces == 0:
        newdata = data
    elif ntraces >= tottraces:
        newdata=np.matrix.mean(data,1) 
    else:
        newdata = np.asmatrix(np.zeros(data.shape))    
        halfwid = int(np.ceil(ntraces/2.0))
        
        # First few traces, that all have the same average
        newdata[:,0:halfwid+1] = np.matrix.mean(data[:,0:halfwid+1],1)
        
        # For each trace in the middle
        for tr in tqdm(range(halfwid,tottraces-halfwid+1)):   
            winstart = int(tr - halfwid)
            winend = int(tr + halfwid)
            newdata[:,tr] = np.matrix.mean(data[:,winstart:winend+1],1) 

        # Last few traces again have the same average    
        newdata[:,tottraces-halfwid:tottraces+1] = np.matrix.mean(data[:,tottraces-halfwid:tottraces+1],1)

    print('done with profile smoothing')
    return newdata, newProfilePos



def tpowGain(data,twtt,power):
    '''
    Apply a t-power gain to each trace with the given exponent.

    INPUT:
    data      data matrix whose columns contain the traces
    twtt      two-way travel time values for the rows in data
    power     exponent

    OUTPUT:
    newdata   data matrix after t-power gain
    '''
    factor = np.reshape(twtt**(float(power)),(len(twtt),1))
    factmat = matlib.repmat(factor,1,data.shape[1])  
    return np.multiply(data,factmat)



def agcGain(data,window):
    '''
    Apply automated gain controll (AGC) by normalizing the energy
    of the signal over a given window width in each trace

    INPUT:
    data       data matrix whose columns contain the traces
    window     window width [in "number of samples"]
    
    OUTPUT:
    newdata    data matrix after AGC gain
    '''
    
    eps=1e-8
    totsamps = data.shape[0]
    # If window is a ridiculous value
    if (window>totsamps):
        # np.maximum is exactly the right thing (not np.amax or np.max)
        energy = np.maximum(np.linalg.norm(data,axis=0),eps)
        # np.divide automatically divides each row of "data"
        # by the elements in "energy"
        newdata = np.divide(data,energy)
    else:
        # Need to go through the samples
        newdata = np.asmatrix(np.zeros(data.shape))
        halfwid = int(np.ceil(window/2.0))
        # For the first few samples, it will always be the same
        energy = np.maximum(np.linalg.norm(data[0:halfwid+1,:],axis=0),eps)
        newdata[0:halfwid+1,:] = np.divide(data[0:halfwid+1,:],energy)
        
        for smp in tqdm(range(halfwid,totsamps-halfwid+1)):
            winstart = int(smp - halfwid)
            winend = int(smp + halfwid)
            energy = np.maximum(np.linalg.norm(data[winstart:winend+1,:],axis=0),eps)
            newdata[smp,:] = np.divide(data[smp,:],energy)

        # For the first few samples, it will always be the same
        energy = np.maximum(np.linalg.norm(data[totsamps-halfwid:totsamps+1,:],axis=0),eps)
        newdata[totsamps-halfwid:totsamps+1,:] = np.divide(data[totsamps-halfwid:totsamps+1,:],energy)          
    return newdata
        

def prepTopo(topofile,delimiter=',',xStart=0):
    '''
    Reads an ASCII text file containing either profile/topo coordinates 
    (if given as two columns) or x,y,z or Easting,Northing,Elevation
    (if given as three columns)

    INPUT:
    topofile    file name for the ASCII text file
    delimiter   delimiter by which the entries are separated 
                (e.g. ',' or tab '\t') [default: ',']
    xStart      if three-dimensional topo data is given:
                profile position of the first x,y,z entry
                [default: 0]

    OUTPUT:
    topoPos     the along-profile coordinates for the elevation points      
    topoVal     the elevation values for the given profile coordinates
    threeD      n x 3 matrix containing the x, y, z values for the 
                topography points
    '''
    
    # Read topofile, see if it is two columns or three columns.
    # Here I'm using numpy's loadtxt. There are more advanced readers around
    # but this one should do for this simple situation
    topotable = np.loadtxt(topofile,delimiter=delimiter)
    topomat = np.asmatrix(topotable)
    # Depending if the table has two or three columns,
    # need to treat it differently
    if topomat.shape[1] is 3:
        # Save the three columns
        threeD = topomat
        # Turn the three-dimensional positions into along-profile
        # distances
        topoVal = topomat[:,2]
        npos = topomat.shape[0]
        steplen = np.sqrt(
            np.power( topomat[1:npos,0]-topomat[0:npos-1,0] ,2.0) + 
            np.power( topomat[1:npos,1]-topomat[0:npos-1,1] ,2.0) +
            np.power( topomat[1:npos,2]-topomat[0:npos-1,2] ,2.0)
        )
        alongdist = np.cumsum(steplen)
        topoPos = np.append(xStart,alongdist+xStart)
    elif topomat.shape[1] is 2:
        threeD = None
        topoPos = topomat[:,0]
        topoVal = topomat[:,1]
        topoPos = np.squeeze(np.asarray(topoPos))
    else:
        print("Something is wrong with the topogrphy file")
        topoPos = None
        topoVal = None
        threeD = None
    return topoPos, topoVal, threeD




def correctTopo(data, velocity, profilePos, topoPos, topoVal, twtt):
    '''
    Corrects for topography along the profile by shifting each 
    Trace up or down depending on provided coordinates.

    INPUT:
    data          data matrix whose columns contain the traces
    velocity      subsurface RMS velocity in m/ns
    profilePos    along-profile coordinates of the traces
    topoPos       along-profile coordinates for provided elevation
                  in meters
    topoVal       elevation values for provided along-profile 
                  coordinates, in meters
    twtt          two-way travel time values for the samples, in ns

    OUTPUT:
    newdata       data matrix with shifted traces, padded with NaN 
    newtwtt       twtt for the shifted / padded data matrix
    maxElev       maximum elevation value
    minElev       minimum elevation value
    '''
    # We assume that the profilePos are the correct along-profile
    # points of the measurements (they can be correted with adj profile)
    # For some along-profile points, we have the elevation from prepTopo
    # So we can just interpolate    
    if not ((all(np.diff(topoPos)>0)) or  (all(np.diff(topoPos)<0))):
        raise ValueError('\x1b[1;31;47m' + 'The profile vs topo file does not have purely increasing or decreasing along-profile positions' + '\x1b[0m')        
    else:
        elev = interp.pchip_interpolate(topoPos,topoVal,profilePos)
        elevdiff = elev-np.min(elev)
        # Turn each elevation point into a two way travel-time shift.
        # It's two-way travel time
        etime = 2*elevdiff/velocity
        timeStep=twtt[3]-twtt[2]
        # Calculate the time shift for each trace
        tshift = (np.round(etime/timeStep)).astype(int)
        maxup = np.max(tshift)
        # We want the highest elevation to be zero time.
        # Need to shift by the greatest amount, where  we are the lowest
        tshift = np.max(tshift) - tshift
        # Make new datamatrix
        newdata = np.empty((data.shape[0]+maxup,data.shape[1]))
        newdata[:] = np.nan
        # Set new twtt
        newtwtt = np.arange(0, twtt[-1] + maxup*timeStep, timeStep)
        nsamples = len(twtt)
        # Enter every trace at the right place into newdata
        for pos in range(0,len(profilePos)):
            #print(type(tshift[pos][0]))
            newdata[tshift[pos][0]:tshift[pos][0]+nsamples ,pos] = np.squeeze(data[:,pos])
        return newdata, newtwtt, np.max(elev), np.min(elev)

    
    

    
def prepVTK(profilePos,gpsmat=None,smooth=True,win_length=51,porder=3):
    '''
    Calculates the three-dimensional coordinates for each trace
    by interpolating the given three dimensional points along the
    profile.

    INPUT:
    profilePos    the along-profile coordinates of the traces
    gpsmat        n x 3 matrix containing the x, y, z coordinates 
                  of given three-dimensional points for the profile
    smooth        Want to smooth the profile's three-dimensional alignment
                  instead of piecewise linear? [Default: True]
    win_length    If smoothing, the window length for 
                  scipy.signal.savgol_filter [default: 51]
    porder        If smoothing, the polynomial order for
                  scipy.signal.savgol_filter [default: 3]

    OUTPUT:
    x, y, z       three-dimensional coordinates for the traces
    '''
    
    if gpsmat is None:
        x = profilePos
        y = np.zeros(x.size)
        z = np.zeros(x.size)
    else:
        #gpstable = np.loadtxt(gpsfile,delimiter=delimiter)
        #gpsmat = np.asmatrix(gpstable)
        #gpsmat=np.asmatrix(gpsmat)
        # Turn the three-dimensional positions into along-profile
        # distances
        if gpsmat.shape[1] is 3:
            npos = gpsmat.shape[0]
            steplen = np.sqrt(
                np.power( gpsmat[1:npos,0]-gpsmat[0:npos-1,0] ,2.0) + 
                np.power( gpsmat[1:npos,1]-gpsmat[0:npos-1,1] ,2.0) +
                np.power( gpsmat[1:npos,2]-gpsmat[0:npos-1,2] ,2.0)
            )
            alongdist = np.cumsum(steplen)
            # gpsPos = np.append(0,alongdist)
            gpsPos = np.append(0,alongdist) + np.min(profilePos)
            # We assume that the profilePos are the correct along-profile
            # points of the measurements (they can be correted with adj profile)
            # For some along-profile points, we have the elevation from prepTopo
            # So we can just interpolate
            xval = gpsmat[:,0]
            yval = gpsmat[:,1]
            zval = gpsmat[:,2]                        
            x = interp.pchip_interpolate(gpsPos,xval,profilePos)
            y = interp.pchip_interpolate(gpsPos,yval,profilePos)
            z = interp.pchip_interpolate(gpsPos,zval,profilePos)
        else:
            npos = gpsmat.shape[0]
            steplen = np.sqrt(
                np.power( gpsmat[1:npos,0]-gpsmat[0:npos-1,0] ,2.0) + 
                np.power( gpsmat[1:npos,1]-gpsmat[0:npos-1,1] ,2.0)  
            )
            alongdist = np.cumsum(steplen)
            # gpsPos = np.append(0,alongdist)
            gpsPos = np.append(0,alongdist) + np.min(profilePos)
            xval = gpsmat[:,0]
            zval = gpsmat[:,1]
            x = interp.pchip_interpolate(gpsPos,xval,profilePos)
            z = interp.pchip_interpolate(gpsPos,zval,profilePos)
            y = np.zeros(len(x))
            
        # Do some smoothing
        if smooth:
            win_length = min(int(len(x)/2),win_length)
            porder = min(int(np.sqrt(len(x))),porder)
            x = signal.savgol_filter(x.squeeze(), window_length=win_length,
                                     polyorder=porder)
            y = signal.savgol_filter(y.squeeze(), window_length=win_length,
                                     polyorder=porder)
            z = signal.savgol_filter(z.squeeze(), window_length=win_length,
                                     polyorder=porder) 
    return x,y,z


def linStackedAmplitude(data,profilePos,twtt,vVals,tVals,typefact):
    '''
    Calculates the linear stacked amplitudes for each two-way 
    travel time sample and the provided velocity range 
    by summing the pixels of the data that follow a line given 
    by the two-way travel time zero offset and the velocity.

    INPUT:
    data          data matrix whose columns contain the traces
    profilePos    along-profile coordinates of the traces
    twtt          two-way travel time values for the samples, in ns
    vVals         list of velocity values for which to calculate the
                  linear stacked amplitudes, in m/ns
    tVals         list of twtt zero-offsets for which to calculate
                  the linear stacked amplitudes, in ns
    typefact      factor for antenna separation depending if this is
                  for CMP (typefact=2) or WARR (typefact=1) data

    OUTPUT:
    linStAmp      matrix containing the linear stacked amplitudes
                  for the given data, tVals, and vVals
    '''
    linStAmp=np.zeros((len(tVals),len(vVals)))
    for vi in tqdm(range(0,len(vVals))):       
        for ti in range(0,len(tVals)):
            t = tVals[ti] + typefact*profilePos/vVals[vi]
            tindices = (np.round((t-twtt[0])/(twtt[3]-twtt[2]))).astype(int)
            # The tindices will be sorted, can use searchsorted because
            # the wave doesn't turn around
            maxi = np.searchsorted(tindices,len(twtt))
            pixels = data[(tindices[0:maxi],np.arange(0,maxi))]
            linStAmp[ti,vi]=np.abs(np.sum(pixels)/pixels.shape[1])
    return linStAmp


def hypStackedAmplitude(data,profilePos,twtt,vVals,tVals,typefact):
    '''
    Calculates the hyperbolic stacked amplitudes for each two-way 
    travel time sample and the provided velocity range 
    by summing the pixels of the data that follow a hyperbola given 
    by the two-way travel time apex and the velocity.

    INPUT:
    data          data matrix whose columns contain the traces
    profilePos    along-profile coordinates of the traces
    twtt          two-way travel time values for the samples, in ns
    vVals         list of velocity values for which to calculate the
                  hyperbolic stacked amplitudes, in m/ns
    tVals         list of twtt zero-offsets for which to calculate
                  the hyperbolic stacked amplitudes, in ns
    typefact      factor for antenna separation depending if this is
                  for CMP (typefact=2) or WARR (typefact=1) data

    OUTPUT:
    hypStAmp      matrix containing the hyperbolic stacked amplitudes
                  for the given data, tVals, and vVals
    '''
    hypStAmp=np.zeros((len(tVals),len(vVals)))
    x2 = np.power(typefact*profilePos,2.0)
    for vi in tqdm(range(0,len(vVals))):       
        for ti in range(0,len(tVals)):
            t = np.sqrt(x2 + 4*np.power(tVals[ti]/2.0 * vVals[vi],2.0))/vVals[vi]
            tindices = (np.round((t-twtt[0])/(twtt[3]-twtt[2]))).astype(int)
            # The tindices will be sorted, can use searchsorted because
            # the wave doesn't turn around
            maxi = np.searchsorted(tindices,len(twtt))
            pixels = data[(tindices[0:maxi],np.arange(0,maxi))]
            hypStAmp[ti,vi]=np.abs(np.sum(pixels)/pixels.shape[1])
    return hypStAmp





# ##### Some helper functions
# def nextpow2(i):
#     n = 1
#     while n < i: n *= 2
#     return n


# def padMat(mat,nrow,ncol):
#     padheight=nrow-mat.shape[0]
#     padwidth=ncol-mat.shape[1]
#     if padheight>0: 
#         mat = np.concatenate((mat,np.zeros((padheight,mat.shape[1]))))
#     if padwidth>0:    
#         pad = np.zeros((nrow,padwidth))
#         mat = np.concatenate((mat,pad),axis=1)
#     return mat


# def padVec(vec,totlen):
#     padwidth=totlen-len(vec)
#     if padwidth>0:
#         vec = np.append(vec,np.zeros(padwidth))
#     return vec


##### Testing / trying to improve performance:
def linStackedAmplitude_alt1(data,profilePos,twtt,vVals,tVals,typefact):
    '''
    Calculates the linear stacked amplitudes for each two-way 
    travel time sample and the provided velocity range 
    by summing the pixels of the data that follow a line given 
    by the two-way travel time zero offset and the velocity.

    INPUT:
    data          data matrix whose columns contain the traces
    profilePos    along-profile coordinates of the traces
    twtt          two-way travel time values for the samples, in ns
    vVals         list of velocity values for which to calculate the
                  linear stacked amplitudes, in m/ns
    tVals         list of twtt zero-offsets for which to calculate
                  the linear stacked amplitudes, in ns
    typefact      factor for antenna separation depending if this is
                  for CMP (typefact=2) or WARR (typefact=1) data

    OUTPUT:
    linStAmp      matrix containing the linear stacked amplitudes
                  for the given data, tVals, and vVals
    '''
    linStAmp=np.zeros((len(tVals),len(vVals)))
    f = interp.interp2d(profilePos, twtt, data)        
    for vi in  tqdm(range(0,len(vVals))):
        for ti in range(0,len(tVals)):
            t = tVals[ti] + typefact*profilePos/vVals[vi]            
            vals = np.diagonal(np.asmatrix(f(profilePos, t)))
            linStAmp[ti,vi] = np.abs(sum(vals)/len(vals))
    return linStAmp


def linStackedAmplitude_alt2(data,profilePos,twtt,vVals,tVals,typefact):
    '''
    Calculates the linear stacked amplitudes for each two-way 
    travel time sample and the provided velocity range 
    by summing the pixels of the data that follow a line given 
    by the two-way travel time zero offset and the velocity.

    INPUT:
    data          data matrix whose columns contain the traces
    profilePos    along-profile coordinates of the traces
    twtt          two-way travel time values for the samples, in ns
    vVals         list of velocity values for which to calculate the
                  linear stacked amplitudes, in m/ns
    tVals         list of twtt zero-offsets for which to calculate
                  the linear stacked amplitudes, in ns
    typefact      factor for antenna separation depending if this is
                  for CMP (typefact=2) or WARR (typefact=1) data

    OUTPUT:
    linStAmp      matrix containing the linear stacked amplitudes
                  for the given data, tVals, and vVals
    '''    
    linStAmp=np.zeros((len(tVals),len(vVals)))
    
    tVals = np.asmatrix(tVals).transpose()   
    for vi in tqdm(range(0,len(vVals))):
        t = tVals + typefact*profilePos/vVals[vi]
        tindices = (np.round((t-twtt[0])/(twtt[3]-twtt[2]))).astype(int)
        for ti in range(0,len(tVals)):
            # The tindices will be sorted, can use searchsorted because
            # the wave doesn't turn around           
            maxi = np.searchsorted(np.ravel(tindices[ti,:]),len(twtt))
            pixels = data[(tindices[ti,0:maxi],np.arange(0,maxi))]
            linStAmp[ti,vi]=np.abs(np.sum(pixels)/pixels.shape[1])
    return linStAmp
