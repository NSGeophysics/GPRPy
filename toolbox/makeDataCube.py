import gprpy as gp
import numpy as np


def reduceSampling(gpr,nprofile,ntwtt):

    samplewidth = int(np.round(gpr.data.shape[1]/nprofile))
    dx = gpr.profilePos[1]-gpr.profilePos[0]
    # First reduce along profile
    datared = np.asarray(np.zeros((gpr.twtt.shape[0]+1,nprofile)))
    profilePosred = np.asarray(np.zeros(nprofile))
    for i in range(0,nprofile):
        datared[:,i] = np.sum(gpr.data[:,i*samplewidth:(i+1)*samplewidth],1)/samplewidth
        profilePosred[i]=np.mean(gpr.profilePos[i*samplewidth:(i+1)*samplewidth])
    gpr.data = datared
    print(gpr.profilePos)
    print(profilePosred)
    
    gpr.profilePos = profilePosred

    ##### To do: Make sure that the averages are over the correct bands and that the
    #####        profile positions are correct too (in the middle of the bands).
    
    
    # Now reduce along twtt
    datared = np.asarray(np.zeros((ntwtt,nprofile)))
    twttred = np.asarray(np.zeros(ntwtt))
    for i in range(0,ntwtt):
        pass
    
    return gpr


def makeDataCube(datalist,outname,nx=50,ny=50,nz=50,nprofile=50,ndepth=50):
    # nprof, ndepth: reduce along profile and time
        
    
    npoints=0
    # Read in all the data points and their topos
    for i in range(0,len(datalist)):
        # These need to have a topo correction
        gpr=gp.gprpy2d(datalist[i])

        gpr=reduceSampling(gpr,nprofile,ndepth)
        
        x,y,z = tools.prepVTK(gpr.profilePos,gpr.threeD)
        X = np.asarray([x.squeeze()])
        Y = np.asarray([y.squeeze()])
        Z = np.reshape(z,(len(z),1)) - np.reshape(gpr.depth,(1,len(gpr.depth)))
        ZZ = np.tile(np.reshape(Z,(1,Z.shape[0],Z.shape[1])),(1,1,1))
        XX = np.tile(np.reshape(X,(X.shape[0],X.shape[1],1)),(1,1,ZZ.shape[2]))
        YY = np.tile(np.reshape(Y,(Y.shape[0],Y.shape[1],1)),(1,1,ZZ.shape[2]))
        data = np.asarray(gpr.data_pretopo.transpose())
        data = np.reshape(data,(1,data.shape[0],data.shape[1]))
        datapoints[i] = np.asarray([np.asarray(XX).flatten(),
                                    np.asarray(YY).flatten(),
                                    np.asarray(ZZ).flatten()]).transpose()
        npoints = npoints+datapoints.shape[1]

    print(npoints)

    
