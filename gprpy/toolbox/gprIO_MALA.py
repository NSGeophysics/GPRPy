import struct
import numpy as np
import re # Regular expressions


def readMALA(file_name):
    # First read header
    info = readGPRhdr(file_name+'.rad')
    try:
        filename = file_name + '.rd3'
        data = np.fromfile(filename, dtype=np.int16)        
    except:
        # I'm not sure what the format of rd7 is. Just assuming it's the same
        filename = file_name + '.rd7'
        data = np.fromfile(filename, dtype=np.int16)
    
    nrows=int(len(data)/int(info['SAMPLES']))
    
    data = (np.asmatrix(data.reshape(nrows,int(info['SAMPLES'])))).transpose()
        
    return data, info
    


def readGPRhdr(filename):
    # Read in text file
    info = {}
    with open(filename) as f:
        for line in f:
            strsp = line.split(':')
            info[strsp[0]] = strsp[1].rstrip()
    return info

