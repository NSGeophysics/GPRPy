import struct
import numpy as np
import array

def readdt1(filename):
    # This function is a python translation of dt1read.m from
    # http://www.lucabaradello.it/files/dt1read.m
    headlen = 32
    with open(filename,"rb") as datafile:
        datafile.seek(8,0) # 0 is beginning of file
        samples, = struct.unpack('f',datafile.read(4))
        samples = int(samples)
        dimtrace = samples*2+128
        datafile.seek(-dimtrace,2) # 2 stands for end of file
        #print(datafile.tell())
        max_traces, = struct.unpack('f',datafile.read(4))
        max_traces = int(max_traces)
        # Initialize matrix
        data = np.zeros((samples,max_traces))
        head = np.zeros((headlen,max_traces))
        # Set the reader to the beginning of the file
        datafile.seek(0,0)
        for j in range(0,max_traces):
            # Read the header info
            for k in range(0,headlen):
                info, =  struct.unpack('f',datafile.read(4))
                head[k,j] = info
            # Now the actual data
            for k in range(0,samples):
                # 2 is the size of short, 'h' is the symbol
                pnt, = struct.unpack('h',datafile.read(2))
                data[k,j] = pnt
            datafile.seek(dimtrace*(j+1),0) 

    return data

        
       
        
