import struct
import numpy as np
import os

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
    # First read header, normalize the file name (get rid of extension)
    file_name, extension = os.path.splitext(file_name)

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
            info[strsp[0]] = strsp[1].strip()

    # if the distance interval is zero, set it to one
    # TODO should be done properly with the coordinates (if available)
    # Alain: I changed 0.1 to eps, in case someone uses high spatial resolution (e.g. lab)
    if float(info['DISTANCE INTERVAL']) < np.finfo(float).eps:
        info['DISTANCE INTERVAL'] = 1.

    return info

