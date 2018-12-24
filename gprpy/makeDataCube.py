import gprpy.gprpy as gp
import gprpy.toolbox.gprpyTools as tools
import numpy as np
import copy
import scipy.interpolate as interp
from scipy.interpolate import griddata
from pyevtk.hl import gridToVTK
from tqdm import tqdm
from scipy.ndimage import gaussian_filter


def reduceSampling(gpr,nprofile,ntwtt):
    
    #gpr2 = copy.copy(gpr)
    if gpr.data_pretopo is None:
        data = gpr.data
        twtt = gpr.twtt       
    else:
        data = gpr.data_pretopo
        twtt = gpr.twtt_pretopo

    if nprofile is None:
        nprofile = data.shape[1]

    if ntwtt is None:
        ntwtt = twtt.shape[0]
        
    samplewidth = int(np.round(data.shape[1]/nprofile))
    nprofile = int(np.ceil(data.shape[1]/samplewidth))

    # This could also be done using scipy's zoom function
    # Both, aplong profile and twtt at the same time
    # First reduce along profile
    datared = np.asarray(np.zeros((twtt.shape[0],nprofile)))
    profilePosred = np.asarray(np.zeros(nprofile))
    for i in range(0,nprofile):
        datared[:,i] = np.mean(data[:,i*samplewidth:(i+1)*samplewidth],1).flatten()
        profilePosred[i]=np.mean(gpr.profilePos[i*samplewidth:(i+1)*samplewidth])
    data = datared    
    gpr.profilePos = profilePosred
    
    # Now reduce along twtt
    samplewidth = int(np.round(data.shape[0]/ntwtt))
    ntwtt = int(np.ceil(data.shape[0]/samplewidth))
    datared = np.asarray(np.zeros((ntwtt,nprofile)))
    twttred = np.asarray(np.zeros(ntwtt))
    for i in range(0,ntwtt):
        datared[i,:] = np.mean(data[i*samplewidth:(i+1)*samplewidth],0)
        twttred[i] = np.mean(twtt[i*samplewidth:(i+1)*samplewidth])
    if gpr.data_pretopo is None:
        gpr.data = datared
        gpr.twtt = twttred
        gpr.depth = gpr.twtt*gpr.velocity/2.0
    else:
        gpr.data_pretopo = datared
        gpr.twtt_pretopo = twttred
        gpr.depth = twttred*gpr.velocity/2.0
    
    return gpr,nprofile,ntwtt



def makeDataCube(datalist,outname,nx=50,ny=50,nz=50,smooth=1,nprofile=None,ndepth=None,method='nearest'):
    # nprofile, ndepth: reduce along profile and time
    
    gpr=gp.gprpy2d(datalist[0])

    gpr,nprofile,ndepth = reduceSampling(gpr,nprofile,ndepth)
            
    # Allocate memory based on nprofile and ndepth. May be overallocating
    allpoints = np.zeros((nprofile*ndepth*len(datalist),3))
    alldata = np.zeros(nprofile*ndepth*len(datalist))
    datalength = np.zeros(len(datalist),dtype=int)
    
    topopoints = 2*np.zeros((nprofile*len(datalist),3))
    topolength = np.zeros(len(datalist),dtype=int)
    
    npoints=0
    # Read in all the data points and their topos
    print('Reading in profiles ...')
    for i in tqdm(range(0,len(datalist))):
        # These need to have a topo correction       
        gpr=gp.gprpy2d(datalist[i])                
        gpr,nprofile,ndepth=reduceSampling(gpr,nprofile,ndepth)
        if i==0:
            currentmaxdepth = np.max(np.abs(gpr.depth))
            depth = gpr.depth
        
        if gpr.data_pretopo is None:
            datalength[i] = gpr.data.shape[0]*gpr.data.shape[1]
        else:
            datalength[i] = gpr.data_pretopo.shape[0]*gpr.data_pretopo.shape[1]
            
        x,y,z = tools.prepVTK(gpr.profilePos,gpr.threeD,smooth=False)
        topolength[i] = len(x)        
        
        Z = np.reshape(z,(len(z),1)) - np.reshape(gpr.depth,(1,len(gpr.depth)))
        
        if np.max(np.abs(gpr.depth)) < currentmaxdepth:
            depth = gpr.depth
            currentmaxdepth = np.max(np.abs(gpr.depth))

        X = np.tile(x,Z.shape[1])
        Y = np.tile(y,Z.shape[1])
        
        indices = np.asarray(np.arange(np.sum(datalength[0:i]),np.sum(datalength[0:i+1])))

        topoindices = np.asarray(np.arange(np.sum(topolength[0:i]),np.sum(topolength[0:i+1])))

        allpoints[indices,:] = np.asarray([X.flatten(),
                                           Y.flatten(),
                                           Z.flatten()]).transpose()
        
        topopoints[topoindices,:] = np.asarray([x,y,z]).squeeze().transpose()
            
        if gpr.data_pretopo is None:
            data = np.asarray(gpr.data.transpose())
            #data = np.asarray(gpr.data)
        else:
            data = np.asarray(gpr.data_pretopo.transpose())
            #data = np.asarray(gpr.data_pretopo)
        
        alldata[indices] = np.reshape(data,(data.shape[0]*data.shape[1]))
        
    # Remove overallocation
    allpoints = allpoints[0:np.sum(datalength),:]
    alldata = alldata[0:np.sum(datalength)]
    topopoints = topopoints[0:np.sum(topolength),:]
    
    # Interpolate
    xg = np.linspace(np.min(allpoints[:,0]),np.max(allpoints[:,0]),nx)
    yg = np.linspace(np.min(allpoints[:,1]),np.max(allpoints[:,1]),ny)
    dg = np.linspace(np.min(depth),np.max(depth),nz)
    [Xg,Yg] = np.meshgrid(xg,yg)
    
    topo = interp.griddata(topopoints[:,0:2],topopoints[:,2],np.asarray([Xg.flatten(),Yg.flatten()]).transpose(),method=method)
    topo = np.reshape(topo,Xg.shape)

    Zg = np.reshape(topo,(topo.shape[0],topo.shape[1],1)) - np.reshape(dg,(1,1,len(dg)))

    
    XXg = (Xg.reshape((Xg.shape[0],Xg.shape[1],1)))*(np.ones((1,1,len(dg))))
    YYg = (Yg.reshape((Yg.shape[0],Yg.shape[1],1)))*(np.ones((1,1,len(dg))))

    intpoints = np.asarray([XXg.flatten(),
                            YYg.flatten(),
                            Zg.flatten()]).transpose()

    print('Interpolating data')
    dataG = interp.griddata(allpoints,alldata,
                            intpoints,
                            method=method)
   
    DG = np.reshape(dataG,  XXg.shape)

    # Smooth
    if smooth is not None:
        DG = gaussian_filter(DG,smooth)
    
    #gridToVTK(outname,XG,YG,ZG,cellData={'gpr data': DG})
    gridToVTK(outname,XXg,YYg,Zg,pointData={'gpr data': DG})
   
