import struct
import numpy as np
import re # Regular expressions


def readMALA(file_name):
    '''
    Reads the MALA .rd3 data file and the .rad header. Can also be used
    to read .rd7 files but I'm not sure if they are really organized
    the same way.

    INPUT: 
    file_name     data file name without the extension!

    OUTPUT:
    data          data matrix whose columns contain the traces
    info          dict with information from the header
    '''
    # First read header
    info = readGPRhdr(file_name+'.rad')
    try:
        filename = file_name + '.rd3'
        data = np.fromfile(filename, dtype=np.int16)        
    except:
        # I'm not sure what the format of rd7 is. Just assuming it's the same
        filename = file_name + '.rd7'
        data = np.fromfile(filename, dtype=np.int32)
    
    nrows=int(len(data)/int(info['SAMPLES']))
    
    data = (np.asmatrix(data.reshape(nrows,int(info['SAMPLES'])))).transpose()
        
    return data, info
    


def readGPRhdr(filename):
    '''
    Reads the MALA header

    INPUT: 
    filename      file name for header with .rad extension
    
    OUTPUT:
    info          dict with information from the header
    '''
    # Read in text file
    info = {}
    with open(filename) as f:
        for line in f:
            strsp = line.split(':')
            info[strsp[0]] = strsp[1].rstrip()
    return info

