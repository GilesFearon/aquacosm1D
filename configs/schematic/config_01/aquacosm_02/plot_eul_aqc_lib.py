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
    z_eul=data_eul.variables['depth'][1:-1] # cropping end points used in eulerian code as boundary conditions - reproduces the croco rho grid by definition 
    time_eul=data_eul.variables['time'][:]/3600/24 # time in days
    chl_eul_avg = np.mean(chl_eul, axis=1)
    data_eul.close()
    return time_eul,z_eul,chl_eul,chl_eul_avg

