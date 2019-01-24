import scipy.interpolate as interp
import numpy as np
from pyevtk.hl import gridToVTK


def interpSurface(pointfile,outfile,nxgrid=100,nygrid=100,method='spline',delimiter='\t',kx=1,ky=1):
    '''
    Creates a surface interpolating the provided three dimensional points.

    INPUT:
    pointfile      ASCII text file with three columns containing 
                   x, y, z or Easting, Northing, Elevation points               
    outfile        filename for VTK file containing the surface 
                   interpolating the given points
    nxgrid         number of mesh points along x-axis
    nygrid         number of mesh points along y-axis
    method         interpolation method: "spline", "nearest",
                   "linear", or "cubic" [default: spline]
    delimiter      for ASCII text input file: what is the delimiter?
                   [default: '\t'  meaining tab]
    kx, ky         If spline interpolation is used: 
                   Spline polynomial order in x and y direction
    '''
    
    # method could be 'nearest', 'linear', 'cubic', or 'spline'.
    # For spline: Set kx and ky. They are the orders of the polynomial
    
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
    if method == 'spline':
        spl = interp.SmoothBivariateSpline(points[:,0],points[:,1],
                                           values,kx=kx,ky=ky)
        Z = spl.ev(X,Y)
    else:
        Z=interp.griddata(points,values,(X,Y),method=method)
        
    colval=np.zeros(Z.shape)
    
    # Write out vts file
    XX = np.tile(X,(1,1,1))
    YY = np.tile(Y,(1,1,1))
    ZZ = np.tile(Z,(1,1,1))  
    gridToVTK(outfile,XX,YY,ZZ)





    

    
