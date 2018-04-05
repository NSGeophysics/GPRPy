import os
#import numpy as np
import gprIO

class gprpy2d:
    def __init__(self,filename=None,desciption=None):
        
        if filename is not None:
            self.importdata(filename)

        
    def importdata(self,filename):
        file_name, file_ext = os.path.splitext(filename)
        if file_ext==".DT1":
            self.data=gprIO.readdt1(filename)
        elif file_ext==".DZT":
            pass
        else:
            print("Can only read dt1 or dzt files")


            
