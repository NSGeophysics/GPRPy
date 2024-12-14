#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 04:35:09 2024

@author: bflinch
"""

import gprpy as gp

testFile = '/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/BotanicalGardens_Oak/Data/DATA07/LINE79.DT1'

mygpr = gp.gprpyProfile()
mygpr.importdata(testFile)
mygpr.dewow(30)
mygpr.remMeanTrace(100)
mygpr.tpowGain(2)

# mygpr.data = bpData(mygpr.data,lf,hf,dt,1)
# #mygpr.topoCorrect('/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/BCZN/Panola/GPRTrench/Line0Topo.txt',delimiter='\t')
# data = np.asarray(mygpr.data)
# #mygpr.depth