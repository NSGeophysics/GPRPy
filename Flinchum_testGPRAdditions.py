#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 04:35:09 2024

@author: bflinch
"""

import gprpy.gprpy as gp
from gprpy.gprpyGUI import GPRPyApp
from gprpy.gprpyCWGUI import GPRPyCWApp
import tkinter as tk


rightcol=9
figrowsp=19+1

root = tk.Tk()

for col in range(rightcol):
    root.columnconfigure(col, weight=1)
for row in range(figrowsp):    
	root.rowconfigure(row, weight=1)

app = GPRPyApp(root)

root.mainloop()



# testFile = '/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/BotanicalGardens_Oak/Data/DATA07/LINE79.DT1'

# mygpr = gp.gprpyProfile()
# mygpr.importdata(testFile)
# mygpr.dewow(1)
# mygpr.remMeanTrace(100)
# mygpr.bpFilter(400, 1200)
# mygpr.plotFreqSpectrum(False)


# mygpr.data = bpData(mygpr.data,lf,hf,dt,1)
# #mygpr.topoCorrect('/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/BCZN/Panola/GPRTrench/Line0Topo.txt',delimiter='\t')
# data = np.asarray(mygpr.data)
# #mygpr.depth