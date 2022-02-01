import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 

from matplotlib.pyplot import *
ion()


for f in pv.GetSources().values():
    pv.Delete(f)

def singMFI(np,pv,dsa,phi,vf):
    print('/mnt/datadrive/syed/allnew/granularEulerFoam/nonlocal/conHopMu'+str(phi)+'vf'+str(vf).replace(".","p")+'/system/controlDict')
    ofreader = pv.OpenFOAMReader(
        FileName = '/mnt/datadrive/syed/allnew/granularEulerFoam/nonlocal/conHopMu'+str(phi)+'vf'+str(vf).replace(".","p")
        +'/system/controlDict') # most clear

    t = np.array(ofreader.TimestepValues)
    N = t.size
    
    # set time to something other than zero to avoid errors about unset fields
    pv.UpdatePipeline(time = t[1])
    
    # if needed, evaluate the point locations used in the simulation
    ofvtkdata = pv.servermanager.Fetch(ofreader)
    ofdata = dsa.WrapDataObject( ofvtkdata)
    ofpts = np.array(ofdata.Points.Arrays[0])
    ptsmin = ofpts.min(axis=0) # minimum values of the three axes
    ptsmax = ofpts.max(axis=0) # maximum values of the three axes
    
    # apply line filter and extract data without rendering anything
    
    pol1 = pv.PlotOverLine(Input=ofreader) #" High Res..." is default?
    pol2 = pv.PlotOverLine(Input=ofreader)
    
    zhalf = (ptsmax[2] - ptsmin[2])/2+ptsmin[2]
    xhalf = (ptsmax[0] - ptsmin[0])/2+ptsmin[0]
    yhalf = (ptsmax[1] - ptsmin[1])/2+ptsmin[1]
    
    Dx = (ptsmax[0] - ptsmin[0])/100
    Dy = (ptsmax[1] - ptsmin[1])/100
    Dz = (ptsmax[2] - ptsmin[2])/100
    
    mu = phi
    D1 = 0.05
    Lw = 0.6
    D2 = 2*(D1/2+Lw*np.sin(mu*np.pi/180))
    
    yM = Lw*np.cos(mu*np.pi/180)
    
    pol1.Source.Point1 = [xhalf, ptsmin[1], zhalf]
    pol1.Source.Point2 = [xhalf, ptsmin[1]+yM, zhalf]
    
    pol2.Source.Point1 = [D1/2-Dx, ptsmin[1], zhalf]
    pol2.Source.Point2 = [D2/2-Dx, ptsmin[1]+yM, zhalf]
    
    
    # extract the line data to numpy array
    
    # get the points for the line
    #vtkdata = pv.servermanager.Fetch(pol) # this returns a vtk object
    poldata1 = dsa.WrapDataObject( pv.servermanager.Fetch(pol1) )
    polpts1 = np.array(poldata1.Points)
    
    x1 = polpts1[:,0]
    #y1 = polpts1[:,1]
    #z1 = polpts1[:,2]
    #poldata2 = dsa.WrapDataObject( pv.servermanager.Fetch(pol2) )
    #polpts2 = np.array(poldata2.Points)
    #x2 = polpts2[:,0]
    #y2 = polpts2[:,1]
    #z2 = polpts2[:,2]
    
    npts = x1.size
    
    # get the values of the simulation data at those points
    up1 = np.zeros((N,npts,3))
    up2 = np.zeros((N,npts,3))
    for i in range(N):
        pv.UpdatePipeline(time = t[i])
        # unfortunately, objected from "Fetch" is not part of the pipeline and must
        # be (re-)called after every update to the time
        poldata1 = dsa.WrapDataObject( pv.servermanager.Fetch(pol1) )
        poldata2 = dsa.WrapDataObject( pv.servermanager.Fetch(pol2) )
        up1[i] = np.array(poldata1.PointData['U.particles'])
        up2[i] = np.array(poldata2.PointData['U.particles'])
        
        
    umean1 = np.mean(up1[:,:,1])
    umean2 = np.mean(up2[:,:,1])
        
    mfi = umean2/umean1
    return mfi 

semiAng = np.array([15,20,25,30,35,40,45])
vf = np.array([0.1,0.2,0.3,0.5,0.7,1])

MFI = np.zeros([semiAng.size,vf.size])

for ii in range(semiAng.size):
    for jj in range(vf.size):
        MFI[ii,jj] = singMFI(np,pv,dsa,semiAng[ii],vf[jj])

np.savez('MFI_data.npz', MFI = MFI)
