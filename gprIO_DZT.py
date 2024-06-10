import struct
import numpy as np
#import re # Regular expressions
from datetime import datetime
import math
from itertools import takewhile

####from contents.py in readgssi softwware ####

# a dictionary of standard gssi antenna codes and frequencies

ANT = {
    '100MHz': 100,
    '200MHz': 200,
    '270MHz': 270,
    '350MHz': 350,
    '400MHz': 400,
    '500MHz': 500,
    '800MHz': 800,
    '900MHz': 900,
    '1600MHz': 1600,
    '2000MHz': 2000,
    '2300MHz': 2300,
    '2600MHz': 2600,
    '3200': 'adjustable',
    '3200MLF': 'adjustable',
    'gprMa': 'adjustable',      # gprMax support
    'GSSI': 'adjustable',       # support for issue #11
    'CUSTOM': 'adjustable',
    '120D' : 120,
    '400D':400,
    '3207': 100,
    '3207AP': 100,
    '5106': 200,
    '5106A': 200,
    '50300': 300,
    '350': 350,
    '350HS': 350,
    'D400HS': 350,
    '50270': 270,
    '50270S': 270,
    'D50300': 300,
    '5103': 400,
    '5103A': 400,
    '50400': 400,
    '50400S': 400,
    '800': 800,
    'D50800': 800,
    '3101': 900,
    '3101A': 900,
    '51600': 1600,
    '51600S': 1600,
    'SS MINI': 1600,
    '62000': 2000,
    '62000-003': 2000,
    '62300': 2300,
    '62300XT': 2300,
    '52600': 2600,
    '52600S': 2600,
}

############################

def readdzt(filename):
    '''
    Reads a GSSI .DZT data file. 

    INPUT: 
    filename     data file name including .DZT extension

    OUTPUT:
    data          data matrix whose columns contain the traces
    info          dict with information from the header

    Thanks to Ian Nesbitt for pointing out extended headers and
    providing the documentation file.
    '''

    # Documentation file is DZT.File.Format.6-14-16.pdf
    
    info = {}
    
    fid = open(filename,'rb');


    # H is unsigned int 16 (ushort = uint16)
    # h is short (int16)
    # I is unsigned int 32 (uint = uint32)
    # i is int32
    # f is float

    # All of the following information is from DZT.File.Format.6-14-16.pdf
    # Provided by Ian Nesbitt
    
    minheadsize = 1024
    infoareasize = 128
    
    info['name'] = fid.name
    
    rh_tag = struct.unpack('h', fid.read(2))[0]  # Pos 00
    
    # Size of the header
    rh_data = struct.unpack('h', fid.read(2))[0] # Pos 02
    
    # Samples per trace
    rh_nsamp = struct.unpack('h', fid.read(2))[0] # Pos 04
    info["Samples per traces :"] = rh_nsamp

    # Bits per word
    rh_bits = struct.unpack('h', fid.read(2))[0] # Pos 06
    
    # Binary offset
    rh_zero = struct.unpack('h', fid.read(2))[0] # Pos 08

    # Scans per second
    rhf_sps = struct.unpack('f', fid.read(4))[0] # Pos 10
    info["Scans per second : "] = rhf_sps
    
    # Scans per meter
    rhf_spm = struct.unpack('f', fid.read(4))[0] # Pos 14
    info["Scans per meter : "] = rhf_spm

    # Meters per mark
    rhf_mpm = struct.unpack('f', fid.read(4))[0] # Pos 18
    info['Meters per mark : '] = rhf_mpm

    # Startposition [ns]
    rhf_position = struct.unpack('f', fid.read(4))[0] # Pos 22
    info["Startposition (ns) : "] = rhf_position

    # length of trace [ns]
    rhf_range = struct.unpack('f', fid.read(4))[0] # Pos 26
    info["Length of trace (ns) : "] = rhf_range
    
    # frequence d echantillonage [GHs]
    info["Sample frequency (GHz) : "]= rh_nsamp/rhf_range

    # Number of passes
    rh_npass = struct.unpack('h', fid.read(2))[0] # Pos 30
    info["Number of passes : "] =  rh_npass

    # Creation date and time
    rhb_cdt = struct.unpack('f', fid.read(4))[0] # Pos 32
    #info["Creation date and time : "] = rhb_cdt
    
    # Last modified date & time
    rhb_mdt = struct.unpack('f', fid.read(4))[0]  # Pos 36
    #info["Last modified date and time : "] = rhb_mdt
    
    # no idea
    rh_mapOffset = struct.unpack('h', fid.read(2))[0] # Pos 40
    #info['rh_mapOffset'] = rh_mapOffset
    # No idea
    rh_mapSize = struct.unpack('h',fid.read(2))[0] # Pos 42
    #info['rh_mapSize'] = rh_mapSize
    
    # offset to text
    rh_text = struct.unpack('h',fid.read(2))[0] # Pos 44

    # Size of text
    rh_ntext = struct.unpack('h',fid.read(2))[0] # Pos 46

    # offset to processing history
    rh_proc = struct.unpack('h',fid.read(2))[0] # Pos 48

    # size of processing history
    rh_nproc = struct.unpack('h',fid.read(2))[0] # Pos 50

    # number of channels
    rh_nchan = struct.unpack('h',fid.read(2))[0] # Pos 52
    info['Number of channels : ']=rh_nchan
    # ... and more stuff we don't really need
    
    #inspired by readgssi software
    fid.seek(98) # start of antenna section
    rh_ant = fid.read(14).decode('utf-8').split('\x00')[0]
    info['Antenna : ']= rh_ant.rsplit('x')[0]
    
    # info['Antenna frequency (MHz) : '] = ANT[info['Antenna : ']]
    try:
        info['Antenna frequency (MHz) : '] = ANT[info['Antenna : ']]
    
    except KeyError:
        print('Antenna frequency couln t be read')

    fid.close()

    # offset will tell us how long the header is in total
    # There could be a minimal header of 1024 bits
    # (very old files may have had 512 bits)
    if rh_data < minheadsize:
        offset = minheadsize*rh_data
    else:
        offset = minheadsize*rh_nchan   

    # Define data type based on words per bit
    # From documentation: Eight byte and sixteen byte samples are
    # unsigned integers. Thirty-two bit samples are signed integers.
    if rh_bits == 8:
        datatype = 'uint8' # unsigned char
    elif rh_bits == 16:
        datatype = 'uint16' # unsigned int
    elif rh_bits == 32:
        datatype = 'int32'

        
    # Read the entire file
    vec = np.fromfile(filename,dtype=datatype)
        
    headlength = offset/(rh_bits/8)

    # Only use the data, discard the header
    datvec = vec[int(headlength):]

    # Turn unsigned integers into signed integers
    # Only necessary where unsigned
    if rh_bits == 8 or rh_bits == 16:
        datvec = datvec - (2**rh_bits)/2.0

    # reshape into matrix
    data = np.reshape(datvec,[int(len(datvec)/rh_nsamp),rh_nsamp])

    data = np.asmatrix(data)
    data_transposed = data.transpose()
    data_transposed[0,:]=0
    data_transposed[1,:]=0
    print('data',data_transposed)
    #print('min, max', np.min(data_transposed),np.max(data_transposed))


    return data_transposed, info