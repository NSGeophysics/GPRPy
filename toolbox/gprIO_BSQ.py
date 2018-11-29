import struct
import numpy as np
import re # Regular expressions


def readBSQ(file_name):
    # First read header
    info = readGPRhdr(file_name+'.GPRhdr')
    filename = file_name + '.dat'
    # Set data type
    if info['data'] == 'float32':
        data = np.fromfile(filename, dtype=np.float32)
    elif info['data'] == 'int16':
        data = np.fromfile(filename, dtype=np.int16)

    data = np.asmatrix(data.reshape( (int(info['lines']), int(info['columns'])) ))
        
    return data, info
    


def readGPRhdr(filename):
    # Read in text file
    info = {}
    with open(filename) as f:
        for line in f:
            strsp = line.split()
            info[strsp[0]] = strsp[-1]
    return info

