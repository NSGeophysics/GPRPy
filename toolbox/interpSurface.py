import scipy.interpolate as interp
import numpy as np
from pyevtk.hl import gridToVTK


def interpSurface(pointfile,outfile,nxgrid=100,nygrid=100,method='nearest',delimiter='\t'):
    # method could be 'nearest' or 'cubic'. To extrapolate: use 'nearest'
    
    # Read in the point values
    pointdata=np.loadtxt(pointfile,delimiter=delimiter)
    points=pointdata[:,0:2]
    values=pointdata[:,2]
    
    # Make interpolation grid
    xmin=np.min(points[:,0])
    xmax=np.max(points[:,0])
    ymin=np.min(points[:,1])
    ymax=np.max(points[:,1])

    x=np.linspace(xmin,xmax,nxgrid)
    y=np.linspace(ymin,ymax,nygrid)
    X,Y=np.meshgrid(x,y)        
    
    # Interpolate the point values
    Z=interp.griddata(points,values,(X,Y),method=method)#'nearest')#'cubic')
    colval=np.zeros(Z.shape)
    
    # Write out vts file
    XX = np.tile(X,(1,1,1))
    YY = np.tile(Y,(1,1,1))
    ZZ = np.tile(Z,(1,1,1))  
    gridToVTK(outfile,XX,YY,ZZ)#,cellData={'surface': colval})





    

    
