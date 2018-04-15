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
    # For each trace
    for tr in tqdm(range(0,data.shape[1])):
        trace = data[:,tr]
        averages = np.zeros(trace.shape)
        # Calculate and subtract running mean
        for i in range(0,data.shape[0]):
            winstart = i - np.floor(window/2.0)
            winend = i + np.floor(window/2.0)
            # If running mean window goes outside of range,
            # set range to "beginning until length"
            if winstart < 0:
                winstart = 0
                winend = window
            # Or to "end-length to end"
            if winend > len(trace):
                winstart = len(trace) - window
                winend = len(trace)
            #print("window for trace %d is %d" %(tr,len(trace[winstart:winend])) )
            #print(np.mean(trace[winstart:winend]))
            #averages[i] = np.mean( trace[winstart:winend])       
            newdata[i,tr] = trace[i] - np.mean(trace[winstart:winend])
    print('done with dewow')
    return newdata

#def remMeanTrace(data):
    

