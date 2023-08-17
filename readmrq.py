#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:34:37 2020

@author: aline
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def read_mrk(filename, nb_trace):
    """
    @author : aline

    Parameters
    ----------
    filename : TYPE string
        DESCRIPTION. path+filename du fichier marqueur à lire
    delta_x : TYPE float ou int
        DESCRIPTION. pas en x (m)
    nb_trace : TYPE int
        DESCRIPTION. nombre de traces

    Returns 
    -------
        TYPE list
        DESCRIPTION. positions (en m) pour chaque trace

    """

    
    fid = open(filename, 'r')
    position = {}
    xp = [] #liste des numeros de trace
    yp = [] #liste des positions (en m)
    
    for line in fid :
        num_trace = line.split()[0]
        pos = line.split()[1]
        
        xp.append(int(num_trace))
        yp.append(int(pos))
    
    fid.close()
    
    delta_x = (yp[-1]-yp[0])/(xp[-1]-xp[0])

    
    f = interpolate.interp1d(xp, yp, fill_value="extrapolate")
    traces = np.linspace(0,nb_trace,nb_trace)
    position = f(traces)
    #plt.plot(position,traces)
    #print('position',position[:nb_trace])

    return position[:nb_trace]