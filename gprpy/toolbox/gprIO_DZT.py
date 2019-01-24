import struct
import numpy as np
#import re # Regular expressions

def readdzt(filename):
    '''
    Reads a GSSI .DZT data file. I learned from A. Tzanis' MatGPR 
    dzt import function how the dzt file header is stored.

    INPUT: 
    filename     data file name including .DZT extension

    OUTPUT:
    data          data matrix whose columns contain the traces
    info          dict with information from the header
    '''
    
    # With help from reviewdztheader.m by Andreas Tzanis

    info = {}
    
    fid = open(filename,'rb');


    # H is unsigned int 16 (ushort = uint16)
    # h is short (int16)
    # I is unsigned int 32 (uint = uint32)
    # i is int32
    # f is float
    
    nchannels = struct.unpack('H', fid.read(2))[0]

    # Size of the header
    headsize = struct.unpack('H', fid.read(2))[0]

    # Samples per trace
    sptrace = struct.unpack('H', fid.read(2))[0]
    info["sptrace"] = sptrace

    # Bits per word
    bpdatum = struct.unpack('H', fid.read(2))[0]

    # Binary offset
    binoffs = struct.unpack('h', fid.read(2))[0]

    # Scans per second
    scpsec = struct.unpack('f', fid.read(4))[0]
    info["scpsec"] = scpsec
    
    # Scans per meter
    scpmeter = struct.unpack('f', fid.read(4))[0]
    info["scpmeter"] = scpmeter

    # Meters per mark
    mpmark = struct.unpack('f', fid.read(4))[0]

    # Startposition
    startposition = struct.unpack('f', fid.read(4))[0]
    info["startposition"] = startposition

    nanosecptrace = struct.unpack('f', fid.read(4))[0]
    info["nanosecptrace"] = nanosecptrace

    scansppass = struct.unpack('H', fid.read(2))[0]
    info["scansppass"] =  scansppass

    fid.close()

    if bpdatum == 8:
        datatype = 'B' # unsigned char
    elif bpdatum == 16:
        datatype = 'H' # unsigned int
    elif bpdatum == 32:
        datatype = 'I'

    # Now read the entire file
    vec=np.fromfile(filename,dtype=datatype)    
        
    # Separate between header and data
    headlength=headsize/(bpdatum/8)
    
    datvec=vec[int(headlength):]

    # Turn unsigned integers into signed integers
    datvec = datvec - (2**bpdatum)/2.

    # reshape into matrix
    #print(sptrace)
    #print(int(len(datvec)/sptrace))
    #data = np.reshape(datvec,[sptrace,int(len(datvec)/sptrace)])
    data = np.reshape(datvec,[int(len(datvec)/sptrace),sptrace])

    data = np.asmatrix(data)
    
    return data.transpose(), info
