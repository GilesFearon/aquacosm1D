import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
ion()

def get_eul_output(eulfile):
    data_eul=Dataset(eulfile)
    chl_eul=data_eul.variables['chl'][:,1:-1] 
    z_eul=data_eul.variables['depth'][1:-1] # cropping end points used in eulerian code as boundary conditions - reproduces the physics rho grid by definition 
    time_eul=data_eul.variables['time'][:]/3600/24 # time in days
    chl_eul_avg = np.mean(chl_eul, axis=1)
    data_eul.close()
    return time_eul,z_eul,chl_eul,chl_eul_avg

def get_aqc_output(aqcfile):
    data_aqc=Dataset(aqcfile)
    p=data_aqc.couping_parameter_p
    chl_aqc=data_aqc.variables['chl'][:,:]
    z_aqc=data_aqc.variables['depth'][:,:]
    z_rank=data_aqc.variables['depth_rank'][:]
    time_aqc=data_aqc.variables['time'][:]/3600/24 # time in days
    chl_aqc_avg = np.mean(chl_aqc, axis=1)
    data_aqc.close()
    return time_aqc,z_aqc,z_rank,chl_aqc,chl_aqc_avg

def get_physics_input(physicsfile):
    
    # get the physics output
    data_physics=Dataset(physicsfile)
    kappa_physics=data_physics.variables['difvho'][:,:,0,0]
    z=data_physics.variables['deptht'][:]
    zw=data_physics.variables['depthw'][:]
    time_physics=data_physics.variables['time_counter'][:]/3600/24 # time in days
    s_flux=data_physics.variables['rsntds'][:,0,0]
    data_physics.close()
    
    return time_physics, z, zw, kappa_physics, s_flux