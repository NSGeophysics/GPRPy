import gprpy as gp
import toolbox.gprpyTools as tools
import numpy as np
import copy
import scipy.interpolate as interp
from pyevtk.hl import gridToVTK


def reduceSampling(gpr,nprofile,ntwtt):
    gpr2 = copy.copy(gpr)
    if gpr.data_pretopo is None:
        data = gpr2.data
        twtt = gpr2.twtt       
    else:
        data = gpr2.data_pretopo
        twtt = gpr2.twtt_pretopo
        
    samplewidth = int(np.round(data.shape[1]/nprofile))
    nprofile = int(np.ceil(data.shape[1]/samplewidth))

    # First reduce along profile
    datared = np.asarray(np.zeros((twtt.shape[0],nprofile)))
    profilePosred = np.asarray(np.zeros(nprofile))
    for i in range(0,nprofile):
        datared[:,i] = np.mean(data[:,i*samplewidth:(i+1)*samplewidth],1).flatten()
        profilePosred[i]=np.mean(gpr2.profilePos[i*samplewidth:(i+1)*samplewidth])
    data = datared    
    gpr2.profilePos = profilePosred
    
    # Now reduce along twtt
    samplewidth = int(np.round(data.shape[0]/ntwtt))
    ntwtt = int(np.ceil(data.shape[0]/samplewidth))
    datared = np.asarray(np.zeros((ntwtt,nprofile)))
    twttred = np.asarray(np.zeros(ntwtt))
    for i in range(0,ntwtt):
        datared[i,:] = np.mean(data[i*samplewidth:(i+1)*samplewidth],0)
        twttred[i] = np.mean(twtt[i*samplewidth:(i+1)*samplewidth])
    if gpr.data_pretopo is None:
        gpr2.data = datared
        gpr2.twtt = twttred
        gpr2.depth = gpr2.twtt*gpr2.velocity/2.0
    else:
        gpr2.data_pretopo = datared
        gpr2.twtt_pretopo = twttred
        gpr2.depth = twttred*gpr2.velocity/2.0
    
    return gpr2


def makeDataCube(datalist,outname,nx=50,ny=50,nz=50,nprofile=None,ndepth=None,method='nearest'):
    # nprofile, ndepth: reduce along profile and time

    gpr=gp.gprpy2d(datalist[0])
    if nprofile is None:
        nprofile = len(gpr.profilePos)

    if ndepth is None:
        if gpr.twtt_pretopo is None:
            ndepth = len(gpr.twtt)
        else:
            ndepth = len(gpr.twtt_pretopo)
            
    # Allocate memory based on nprofile and ndepth. May be overallocating
    allpoints = np.zeros((nprofile*ndepth*len(datalist),3))
    alldata = np.zeros(nprofile*ndepth*len(datalist))
    datalength = np.zeros(len(datalist)+1,dtype=int)
    
    npoints=0
    # Read in all the data points and their topos
    for i in range(0,len(datalist)):
        # These need to have a topo correction       
        gpr=gp.gprpy2d(datalist[i])
        
        gpr=reduceSampling(gpr,nprofile,ndepth)
        if gpr.data_pretopo is None:
            datalength[i+1] = gpr.data.shape[0]*gpr.data.shape[1]
        else:
            datalength[i+1] = gpr.data_pretopo.shape[0]*gpr.data_pretopo.shape[1]

        x,y,z = tools.prepVTK(gpr.profilePos,gpr.threeD,smooth=False)
        
        X = np.asarray([x.squeeze()])
        Y = np.asarray([y.squeeze()])
        Z = np.reshape(z,(len(z),1)) - np.reshape(gpr.depth,(1,len(gpr.depth)))

        XX = np.tile(X,Z.shape[1]).flatten()
        YY = np.tile(Y,Z.shape[1]).flatten()
        ZZ = Z.flatten()        
        
        indices = np.array(np.arange(np.sum(datalength[0:i+1]),np.sum(datalength[0:i+2])))
        
     
        allpoints[indices,:] = np.asarray([XX.transpose(),
                                           YY.transpose(),
                                           ZZ.transpose()]).transpose()
        
        if gpr.data_pretopo is None:
            data = np.asarray(gpr.data.transpose())
        else:
            data = np.asarray(gpr.data_pretopo.transpose())

        
        alldata[indices] = np.reshape(data,(data.shape[0]*data.shape[1]))                
        
    # Remove overallocation
    allpoints = allpoints[0:np.sum(datalength),:]
    alldata = alldata[0:np.sum(datalength)]

    # Interpolate
    xg = np.linspace(np.min(allpoints[:,0]),np.max(allpoints[:,0]),nx)
    yg = np.linspace(np.min(allpoints[:,1]),np.max(allpoints[:,1]),ny)
    zg = np.linspace(np.min(allpoints[:,2]),np.max(allpoints[:,2]),nz)
    XG,YG,ZG = np.meshgrid(xg,yg,zg)

    dataG = interp.griddata(allpoints,alldata,
                            (XG.flatten(),YG.flatten(),ZG.flatten()),
                            method=method)
    
    DG = np.reshape(dataG,XG.shape)
   
    #DG=DG[0:-1,0:-1,0:-1]
    
    gridToVTK(outname,XG,YG,ZG,cellData={'gpr data': DG})

