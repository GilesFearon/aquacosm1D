#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 12:32:12 2021

@author: osboxes
"""

import numpy as np
import xarray as xr

ndays = 100
calendar = np.linspace(0,ndays,num=ndays*48+1)*24.*3600. # in seconds.

mlds = [10,20,50,100]
kappas = [0.00001,0.0001,0.001,0.01]

for mld in mlds:
    for kappa in kappas:
        
        # create the z grid (1 m grid)
        z_w = np.linspace(0, mld, num=mld+1, endpoint=True,)
        z_r = 0.5*(z_w[1:]+z_w[0:-1])
        
        Akt_store    = np.zeros((len(calendar),len(z_w))) + kappa
        
        # dirunal surface radiation
        omega = 2*np.pi/86400. # diurnal frequency in seconds
        Qswmax = 800 # diurnal peak surface radiation
        srflx_store = np.cos(omega*calendar - np.pi)  * Qswmax
        srflx_store[srflx_store<0]=0
        
        # write a netcdf in the same format at the NEMO netcdf file which the aquacosm model is hard coded to ingest

        Akt_store=np.expand_dims(np.expand_dims(Akt_store,axis=-1),axis=-1)
        srflx_store=np.expand_dims(np.expand_dims(srflx_store,axis=-1),axis=-1)
        ds = xr.Dataset(data_vars=dict(difvho=(["time_counter","depthw","y","x"],Akt_store),
                                               rsntds=(["time_counter","y","x"],srflx_store),
                                               ),
                                               coords=dict(depthw=(["depthw"],z_w),deptht=(["deptht"],z_r),time_counter=(["time_counter"],calendar),x=(["x"],[0]),y=(["y"],[0])))
        ds.time_counter.attrs['units'] = 'seconds since 0000-1-1 00:00:00' # dummy time as an analytical experiment
        # write the output
        fname_out="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
        ds.to_netcdf(fname_out)