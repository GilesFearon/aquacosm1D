import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
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
            # terrible coding this, but just a quick fix to get around when temp_thermocline is out of range
            # this can happen when the subsurface warms (due to mixing of warm surface water) to over the spcified thermocline depth
            continue
    z_therm[z_therm == 0] = np.NaN
    return z_therm

def get_Cs_aquacosm(time_eul,z_therm_eul,time_aqc,z_aqc,tpas_aqc):
    # get the mean concentration of aquacosm tracer tracer tpas_aqc above a depth of z_therm_eul
    
    Cs=np.zeros(len(time_aqc))
    
    for t in range(0,len(time_aqc)):
        # find the closest eulerian time
        # (aquacosm time will be contained within the range of the croco time as it had to run with this as input)
        t_eul = (np.abs(time_aqc[t] - time_eul)).argmin()
        #
        try:
            Cs[t]=tpas_aqc[t,z_aqc[t,:]<z_therm_eul[t_eul]].mean()
        except:
            # terrible coding this, but just a quick fix to get around when z_therm[t]=NaN
            continue
    
    Cs[Cs == 0] = np.NaN     
    return Cs

def get_Cs_aquacosm_2(time_eul,z_therm_eul,time_aqc,z_aqc,tpas_aqc):
    # get the mean concentration of aquacosm tracer tracer tpas_aqc above a depth of z_therm_eul
    # loosely based on __integrate_chl__ in aquacosm1D_reactions.py
    # not using this as it gives basically the same thing, but more noisy results
    
    Cs=np.zeros(len(time_aqc))
    
    for t in range(0,len(time_aqc)):
        # find the closest eulerian time
        # (aquacosm time will be contained within the range of the croco time as it had to run with this as input)
        t_eul = (np.abs(time_aqc[t] - time_eul)).argmin()
        #
        # subset the aquacosm data to the surface particles only
        tpas_aqc_surf=tpas_aqc[t,z_aqc[t,:]<z_therm_eul[t_eul]]
        z_aqc_surf=z_aqc[t,z_aqc[t,:]<z_therm_eul[t_eul]]
        #
        # integrate the concentrations weighted by the distance between the particles
        I = np.zeros_like(tpas_aqc_surf)
        I[0]  = tpas_aqc_surf[0]*z_aqc_surf[0] #assume uniform tpas from z=0 down to the first particle
        I[1:] = 0.5*(tpas_aqc_surf[1:]+tpas_aqc_surf[:-1]) * (z_aqc_surf[1:]-z_aqc_surf[:-1])
        #
        # divide by the layer thickness to get the mean
        Cs[t]=np.sum(I)/z_therm_eul[t_eul]
            
    Cs[Cs == 0] = np.NaN     
    return Cs

def get_Cs_eulerian(time,z,zw,tpas,z_therm):
    # get the mean concentration of tracer tpas above a depth of z_therm
    
    Cs=np.zeros(len(time))
    
    z_thickness=zw[1:]-zw[0:-1];
    
    for t in range(0,len(time)):
        try:
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
        except:
            # terrible coding this, but just a quick fix to get around when z_therm[t]=NaN
            continue
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
    
def get_eul_output(eulfile):
    data_eul=Dataset(eulfile)
    # Chl_C=np.float(data_eul.Chl_C) # ratio of chlorophyll to carbon [mgChl/mgC]
    chl_eul=data_eul.variables['chl'][:,1:-1] 
    z_eul=data_eul.variables['depth'][1:-1] # cropping end points used in eulerian code as boundary conditions - produces same grid as croco by definition 
    time_eul=data_eul.variables['time'][:]/3600/24 # time in days
    time_eul=time_eul-time_eul[0]
    data_eul.close()
    
    return time_eul,z_eul,chl_eul

def get_aqc_output(aqcfile):
    
    # get the aquacosm output
    data_aqc=Dataset(aqcfile)
    p=data_aqc.couping_parameter_p
    chl_aqc=data_aqc.variables['chl'][:,:]
    z_aqc=data_aqc.variables['depth'][:,:]
    z_rank=data_aqc.variables['depth_rank'][:]
    time_aqc=data_aqc.variables['time'][:]/3600/24 # time in days
    data_aqc.close()
    return time_aqc,z_aqc,z_rank,chl_aqc,p

def w2rho(time_croco,zw,z,kappa_w):
    # interpolate a w grid kappa value onto the rho vertical grid
    kappa_r=np.zeros((len(time_croco),len(z)))
    for t in range(0,len(time_croco)):
        kappa_r[t,:]=interp1d(zw,kappa_w[t,:],kind='linear')([z])
    return kappa_r


def get_aqc_diags(diagfile):
    # get the aquacosm diagnostics output
    data_aqc=Dataset(diagfile)
    time_aqc=data_aqc.variables['time'][:] # time in days
    z_diag=data_aqc.variables['depth_diag'][:]
    z_aqc=data_aqc.variables['depth_aqc'][:,:]
    rank_aqc=data_aqc.variables['rank_aqc'][:]
    r=data_aqc.variables['r'][:,:]
    stdev_r=data_aqc.variables['stdev_r'][:,:]
    stdev_r_norm=data_aqc.variables['stdev_r_norm'][:,:]
    stdev_chl=data_aqc.variables['stdev_chl'][:,:]
    stdev_chl_norm=data_aqc.variables['stdev_chl_norm'][:,:]
    data_aqc.close()
    return time_aqc,z_diag,z_aqc,rank_aqc,r,stdev_r,stdev_r_norm,stdev_chl,stdev_chl_norm

def get_aqc_reactions(time,z_aqc,z_rank,chl_aqc,React):
    # get the reaction rates normalised by the concentration of the active tracer
        
    # take the chl/C conversion into account, as input C represents chlorophyl
    # but BioShading_onlyC needs Carbon input
    C=chl_aqc/React.Chl_C 
    
    # set up an aquacosm array in the correct format to be given to the reactions library
    Nt,Npts=shape(C) 
    Nscalars=1
    Particles      = Aquacosm1D_Particles(
            zeros((Npts, Nscalars+2), dtype='float64')
            )
    Particles[:,0] = z_rank
    
    # compute the reaction rates at each time-step, normalised by the concentration of the activate tracer
    RRates_normalised=np.zeros_like(C)
    for ii in range(len(time)):
        React.wc.set_current_time(time[ii])
        Particles[:,1]=z_aqc[ii,:]
        Particles[:,2]=C[ii,:]
        RRates=React.current_model(Particles, React.wc, time[ii])
        RRates_normalised[ii,:]=RRates[:,0]/Particles[:,2]
        
    return RRates_normalised


