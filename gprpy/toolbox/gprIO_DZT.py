import struct
import numpy as np
#import re # Regular expressions

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
    
    rh_tag = struct.unpack('h', fid.read(2))[0]  # Pos 00
    
    # Size of the header
    rh_data = struct.unpack('h', fid.read(2))[0] # Pos 02
    
    # Samples per trace
    rh_nsamp = struct.unpack('h', fid.read(2))[0] # Pos 04
    info["rh_nsamp"] = rh_nsamp

    # Bits per word
    rh_bits = struct.unpack('h', fid.read(2))[0] # Pos 06
    
    # Binary offset
    rh_zero = struct.unpack('h', fid.read(2))[0] # Pos 08

    # Scans per second
    rhf_sps = struct.unpack('f', fid.read(4))[0] # Pos 10
    info["rhf_sps"] = rhf_sps
    
    # Scans per meter
    rhf_spm = struct.unpack('f', fid.read(4))[0] # Pos 14
    info["rhf_spm"] = rhf_spm

    # Meters per mark
    rhf_mpm = struct.unpack('f', fid.read(4))[0] # Pos 18

    # Startposition [ns]
    rhf_position = struct.unpack('f', fid.read(4))[0] # Pos 22
    info["rhf_position"] = rhf_position

    # length of trace [ns]
    rhf_range = struct.unpack('f', fid.read(4))[0] # Pos 26
    info["rhf_range"] = rhf_range

    # Number of passes
    rh_npass = struct.unpack('h', fid.read(2))[0] # Pos 30
    info["rh_npass"] =  rh_npass

    # Creation date and time
    rhb_cdt = struct.unpack('f', fid.read(4))[0] # Pos 32
    info["rhb_cdt"] = rhb_cdt
    
    # Last modified date & time
    rhb_mdt = struct.unpack('f', fid.read(4))[0]  # Pos 36
    
    # no idea
    rh_mapOffset = struct.unpack('h', fid.read(2))[0] # Pos 40

    # No idea
    rh_mapSize = struct.unpack('h',fid.read(2))[0] # Pos 42

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

    # ... and more stuff we don't really need

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
    
    return data.transpose(), info
