import numpy as np
import numpy.matlib as matlib
import scipy.interpolate as interp
# For progress bar
import time
from tqdm import tqdm


def timeZeroAdjust(data):
        
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


def remMeanTrace(data,ntraces):
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



def tpowGain(data,twtt,power):
    factor = np.reshape(twtt**(float(power)),(len(twtt),1))
    factmat = matlib.repmat(factor,1,data.shape[1])
    
    return np.multiply(data,factmat)



def agcGain(data,window):
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
        

def prepTopo(topofile,delimiter=','):
    # Read topofile, see if it is two columns or three columns.
    # Here I'm using numpy's loadtxt. There are more advanced readers around
    # but this one should do for this simple situation
    #delimiter = ','
    topotable = np.loadtxt(topofile,delimiter=delimiter)
    topomat = np.asmatrix(topotable)
    # Depending if the table has two or three columns,
    # need to treat it differently
    if topomat.shape[1] is 3:
        # Turn the three-dimensional positions into along-profile
        # distances
        topoVal = topomat[:,2]
        npos = topomat.shape[0]
        steplen = np.sqrt(
            np.power( topomat[1:npos,0]-topomat[0:npos-1,0] ,2.0) + 
            np.power( topomat[1:npos,1]-topomat[0:npos-1,1] ,2.0)
        )
        alongdist = np.cumsum(steplen)
        topoPos = np.append(0,alongdist)
    elif topomat.shape[1] is 2:
        topoPos = topomat[:,0]
        topoVal = topomat[:,1]
        topoPos = np.squeeze(np.asarray(topoPos))
        
    return topoPos, topoVal




def correctTopo(data, velocity, profilePos, topoPos, topoVal, twtt):
    # The variable "topoPos" provides the along-profile coordinates
    # for which the topography is given. 
    # We allow several possibilities to provide topoPos:
    # If None is given, then we assume that they are regularly
    # spaced along the profile
    if topoPos is None:
        topoPos = np.linspace(np.min(profilePos),np.max(profilePos),
                              np.size(topoVal))
    # If it's an integer or a float, then it gives the evenly spaced
    # intervals
    elif type(topoPos) is int:
        topoPos = np.arange(np.min(profilePos),np.max(profilePos),
                            float(topoPos))
    elif type(topoPos) is float:
        topoPos = np.arange(np.min(profilePos),np.max(profilePos),
                            topoPos)
    # Or it could be a file giving the actual along-profile positions    
    elif type(topoPos) is str:
        delimiter = ','
        topopostable = np.loadtxt(topoPos,delimiter=delimiter)

    # Next we need to interpolate the topography
    elev = interp.pchip_interpolate(topoPos,topoVal,profilePos)

    elevdiff = elev-np.min(elev)
    # Turn each elevation point into a two way travel-time shift.
    # It's two-way travel time
    etime = 2*elevdiff/velocity

    timeStep=twtt[1]-twtt[0]
    
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

    return newdata, newtwtt, np.max(elev)

    
    

    
