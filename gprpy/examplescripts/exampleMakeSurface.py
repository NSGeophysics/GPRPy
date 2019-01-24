from gprpy.interpSurface import *

xsamples=200
ysamples=200

# I usually like spline best. It is smooth and extrapolates
interpSurface('../exampledata/pickedSurfaceData/testpick_3D.txt','exampleSurf_skx1ky2',xsamples,ysamples,'spline',kx=1,ky=2)

# 'nearest' extrapolates but makes jumpy surface 
interpSurface('../exampledata/pickedSurfaceData/testpick_3D.txt','exampleSurf_n',xsamples,ysamples,'nearest')

# 'cubic' does not extrapolate but the surface is smooth
interpSurface('../exampledata/pickedSurfaceData/testpick_3D.txt','exampleSurf_c',xsamples,ysamples,'cubic')

# 'linear' is nice too
interpSurface('../exampledata/pickedSurfaceData/testpick_3D.txt','exampleSurf_l',xsamples,ysamples,'linear')

