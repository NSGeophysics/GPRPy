import numpy as np
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
    newdata = np.asmatrix(np.zeros(data.shape))
    
    # If the window is larger or equal to the number of samples,
    # then we can do a much faster dewow
    if (window >= data.shape[0]):
        for tr in tqdm(range(0,data.shape[1])):
            newdata[:,tr]=data[:,tr]-np.mean(data[:,tr])*np.ones((data.shape[0],1))
    else:
        # For each trace
        for tr in tqdm(range(0,data.shape[1])):
            trace = data[:,tr]
            averages = np.zeros(trace.shape)
            # Calculate and subtract running mean
            for i in range(0,data.shape[0]):
                winstart = int(i - np.floor(window/2.0))
                winend = int(i + np.floor(window/2.0))
                # If running mean window goes outside of range,
                # set range to "beginning until length"
                if winstart < 0:
                    winstart = 0
                    winend = window
                # Or to "end-length to end"
                if winend > len(trace):
                    winstart = len(trace) - window
                    winend = len(trace)     
                newdata[i,tr] = trace[i] - np.mean(trace[winstart:winend])
    print('done with dewow')
    return newdata


def remMeanTrace(data,ntraces):
    newdata = np.asmatrix(np.zeros(data.shape))
    tottraces = data.shape[1]

    # For each trace
    for tr in tqdm(range(0,data.shape[1])):
        winstart = int(tr - np.floor(ntraces/2.0))
        winend = int(tr + np.floor(ntraces/2.0))
        if (winstart < 0):
            winstart = 0
            winend = min(ntraces,tottraces)
        elif (winend > tottraces):
            winstart = max(tottraces - ntraces,0)
            winend = tottraces
           
        avgtr = np.zeros(data[:,tr].shape)
        for i in range(winstart,winend):
            avgtr = avgtr + data[:,i]
            
        avgtr = avgtr/float(winend-winstart)
            
        newdata[:,tr] = data[:,tr] - avgtr
            
    return newdata

#def tpowGain(data,twtt,power):
#    factor = np.reshape(twtt**(float(power)),
#                                        (len(twtt),1)))
#    factmat = np.matlib.repmat(factor,(data).shape[1]))
#
#    return np.multiply(data,factmat)
