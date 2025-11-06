import struct
import numpy as np
import re # Regular expressions

def readdt1(filename):
    '''
    Reads the Sensors and Software .DT1 data files. This function is
    a Python translation of http://www.lucabaradello.it/files/dt1read.m

    INPUT: 
    filename      data file name including the .DT1 extension

    OUTPUT:
    data          data matrix whose columns contain the traces
    '''
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
    return np.asmatrix(data)

        
       
def readdt1Header(filename):
    '''
    Reads the Sensors and Software .HD header files.

    INPUT: 
    filename      header file name including the .HD extension

    OUTPUT:
    info          dict with information from the header
    '''
    
    info = {}
    with open(filename,"r",newline='\n') as datafile:
        datafile.readline().strip()
        info["system"] = datafile.readline().strip()
        info["date"] = datafile.readline().strip()
        string = datafile.readline().strip()
        var = re.match(r'NUMBER OF TRACES   = (.*)', string)
        info["N_traces"] = int(var.group(1)) 
        string = datafile.readline().strip()
        var = re.match(r'NUMBER OF PTS/TRC  = (.*)', string)
        info["N_pts_per_trace"] = int(var.group(1)) 
        string = datafile.readline().strip()
        var = re.match(r'TIMEZERO AT POINT  = (.*)', string)
        info["TZ_at_pt"] = float(var.group(1)) 
        string = datafile.readline().strip()
        var = re.match(r'TOTAL TIME WINDOW  = (.*)', string)
        info["Total_time_window"] = float(var.group(1)) 
        string = datafile.readline().strip()
        var = re.match(r'STARTING POSITION  = (.*)', string)
        info["Start_pos"] = float(var.group(1)) 
        string = datafile.readline().strip()
        var = re.match(r'FINAL POSITION     = (.*)', string)
        info["Final_pos"] = float(var.group(1))
        string = datafile.readline().strip()
        var = re.match(r'STEP SIZE USED     = (.*)', string)
        info["Step_size"] = float(var.group(1))
        string = datafile.readline().strip()
        var = re.match(r'POSITION UNITS     = (.*)', string)
        info["Pos_units"] = str(var.group(1))
        string = datafile.readline().strip()
        var = re.match(r'NOMINAL FREQUENCY  = (.*)', string)
        info["Freq"] = float(var.group(1))
        string = datafile.readline().strip()
        var = re.match(r'ANTENNA SEPARATION = (.*)', string)
        info["Antenna_sep"] = float(var.group(1))                   
        # If you need more of the header info, you can just continue as above
        # Transform feet to meters
        if info['Pos_units'] == 'ft':
            info["Start_pos"] = info["Start_pos"]*0.3048
            info["Final_pos"] = info["Final_pos"]*0.3048
            info["Step_size"] = info["Step_size"]*0.3048
            info["Antenna_sep"] = info["Antenna_sep"]*0.3048
            info['Pos_units'] = 'm'
    return info
