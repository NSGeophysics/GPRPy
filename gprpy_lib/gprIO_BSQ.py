import struct
import numpy as np
import re # Regular expressions


def readBSQ(file_name):
    '''
    Reads the ENVI standard BSQ files. The file extension for the 
    data needs to be ".dat" and the extension for the header needs 
    to be ".GPRhdr" 

    INPUT: 
    file_name      data file name without the extension!

    OUTPUT:
    data          data matrix whose columns contain the traces
    info          dict with information from the header
    '''    
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
    '''
    Reads the ENVI standard BSQ file header

    INPUT: 
    filename      file name for header with extension
    
    OUTPUT:
    info          dict with information from the header
    '''
    # Read in text file
    info = {}
    with open(filename) as f:
        for line in f:
            strsp = line.split()
            info[strsp[0]] = strsp[-1]
    return info

