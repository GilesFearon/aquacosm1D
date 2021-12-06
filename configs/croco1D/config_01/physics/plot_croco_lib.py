# import sys
import numpy as np
from netCDF4 import Dataset
# from datetime import datetime, timedelta
# from pathlib import Path
from scipy.interpolate import interp1d
# ion()

def get_z_therm_croco(time,z,temp,temp_thermocline):
    z_therm=np.zeros(len(time))
    for t in range(0,len(time)):
        try:
            z_therm[t]=interp1d(temp[t,:],z,kind='linear')(temp_thermocline)
        except:
            continue
    z_therm[z_therm == 0] = np.NaN
    return z_therm

def get_Cs_eulerian(time,z,zw,tpas,z_therm):
    # get the mean concentration of tracer tpas above a depth of z_therm
    
    Cs=np.zeros(len(time))
    
    z_thickness=zw[1:]-zw[0:-1];
    
    for t in range(0,len(time)):
        
        if not np.isnan(z_therm[t]):
        
            # find the z indices above and below the interpolated depth
            indx_above=np.where(z<z_therm[t])[-1][-1]
            indx_below=indx_above+1
            # compute the weights to apply to the indices either side of the depth of the thermocline
            weight_above=(z[indx_above]-z_therm[t])/(z[indx_above]-z[indx_below]);
            weight_below=1-weight_above;
            #
            # integrate tpas from the surface to z_therm
            for k in range(0,indx_above): # using range means we don't include indx_above in the loop, which is what we want
                Cs[t]=Cs[t] + tpas[t,k]*z_thickness[k]
            # add half of the layer above
            Cs[t]=Cs[t]+0.5*tpas[t,indx_above]*z_thickness[indx_above]
            # add the mean of the layer above and below, weighted by where the
            # thermocline is relative to the two layers
            Cs[t]=Cs[t]+weight_above*np.mean(tpas[t,indx_above:indx_below])*np.mean(z_thickness[indx_above:indx_below])
            # and now divide by the layer thickness to get the mean
            Cs[t]=Cs[t]/z_therm[t]
    
    Cs[Cs == 0] = np.NaN 
       
    return Cs

def get_croco_output(crocofile):
    
    # get the croco output
    data_croco=Dataset(crocofile)
    tpas_croco=data_croco.variables['tpas'][:,:,0,0]
    temp_croco=data_croco.variables['temp'][:,:,0,0]
    kappa_croco=data_croco.variables['difvho'][:,:,0,0]
    u_croco=data_croco.variables['u'][:,:,0,0]
    v_croco=data_croco.variables['v'][:,:,0,0]
    z=data_croco.variables['deptht'][:]
    zw=data_croco.variables['depthw'][:]
    time_croco=data_croco.variables['time_counter'][:]/3600/24 # time in days
    s_flux=data_croco.variables['rsntds'][:,0,0]
    tau_x=data_croco.variables['tau_x'][:,0,0]
    tau_y=data_croco.variables['tau_y'][:,0,0]
    dzetadx=data_croco.variables['dzetadx'][:,0,0]
    data_croco.close()
    
    return time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx

def w2rho(time_croco,zw,z,kappa_w):
    # interpolate a w grid value onto the rho vertical grid
    kappa_r=np.zeros((len(time_croco),len(z)))
    for t in range(0,len(time_croco)):
        kappa_r[t,:]=interp1d(zw,kappa_w[t,:],kind='linear')([z])
    return kappa_r
    
    
