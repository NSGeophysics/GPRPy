import numpy as np



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




#def remMeanTrace(data):
    

