import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
import xarray as xr
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import params # params.py must be in this directory

ion()

def do_the_diags(react_params,physics_params,p):
    
    crocodir='../physics/'
    crocofilename="mean"+str(physics_params.mean_tau)+"_mld"+str(physics_params.mld)+"_amp"+str(physics_params.amplitude)+"_flx"+str(physics_params.Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename
    
    dt = 5             
    wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
    React = set_up_reaction(wc, dt, BioShading_onlyC,
                        LightDecay=react_params.LightDecay,
                        AlphaEpsilon=react_params.AlphaEpsilon,
                        MaxPhotoRate = react_params.MaxPhotoRate, 
                        BasalMetabolism = react_params.BasalMetabolism,
                        Chl_C = react_params.Chl_C,
                        CrowdingMortality = react_params.CrowdingMortality,
                        CrowdingHalfSaturation = react_params.CrowdingHalfSaturation,
                        Chl_light_abs = react_params.Chl_light_abs)
    
    fname='aquacosm_p'+"{0:1.0e}".format(p)+'_'+react_params.Name+'_'+physics_params.Name
    aqcfile=fname+'.nc'
    print('\n working on ' + aqcfile +'\n')
    
    # get the aquacosm data
    time_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
    Nt=len(time_aqc)
    
    # compute the normalised realised reaction
    # directly from the reactions function
    r = get_aqc_reactions(time_aqc*3600*24,z_aqc,rank_aqc,chl_aqc,React)
    r = r*3600*24 # now days^-1
    
    # compute the patchiness of the aquacosms
    # (normalised standard deviation in the vicinity of each aquacosm)
    rad=2.0 # search radius (m) for computation of local stdev around each aquacosm
    
    # create the z grid for diagnostics 
    # using a regular 1 m grid over 50 m (i.e. total water depth)
    z_out = np.linspace(0, 50, num=51, endpoint=True,)
    Np=len(z_out)
    
    # compute the diagnostics
    stdev_chl=np.zeros((Nt,Np))
    stdev_chl_norm=np.zeros((Nt,Np))
    stdev_r=np.zeros((Nt,Np))
    stdev_r_norm=np.zeros((Nt,Np))
    for ii in range(Nt):
        # loop through every diagnostic depth and find particles within search radius to compute stdevs
        for dd in range(Np):
            dist_relative=np.abs(z_out[dd]-z_aqc[ii,:]) # distance of aquacosms from this diagnostic depth
            stdev_chl[ii,dd]=np.std(chl_aqc[ii,dist_relative<rad])
            stdev_chl_norm[ii,dd]=stdev_chl[ii,dd]/np.mean(chl_aqc[ii,dist_relative<rad])
            stdev_r[ii,dd]=np.std(r[ii,dist_relative<rad])
            stdev_r_norm[ii,dd]=stdev_r[ii,dd]/np.mean(r[ii,dist_relative<rad])
    
    # write a netcdf output file
    fname_out=fname+'_diags.nc'
    ds = xr.Dataset(data_vars=dict(depth_aqc=(["time","rank_aqc"],z_aqc),
                                    r=(["time","rank_aqc"],r),
                                    stdev_r=(["time","depth_diag"],stdev_r),
                                    stdev_r_norm=(["time","depth_diag"],stdev_r_norm),
                                    stdev_chl=(["time","depth_diag"],stdev_chl),
                                    stdev_chl_norm=(["time","depth_diag"],stdev_chl_norm),
                                    ),
                                    coords=dict(rank_aqc=(["rank_aqc"],rank_aqc),
                                                depth_diag=(["depth_diag"],z_out),
                                                time=(["time"],time_aqc)))
    # coords metadata
    ds.time.attrs['units'] = 'seconds since start of simulation' # not a real time period as an analytical experiment
    ds.rank_aqc.attrs['long_name'] = 'depth rank of aquacosms'
    ds.rank_aqc.attrs['units'] = '-'
    ds.depth_diag.attrs['long_name'] = 'depth of diagnostic output'
    ds.depth_diag.attrs['units'] = 'm'
    # variable metadata
    ds.depth_aqc.attrs['long_name'] = 'depth of aquacosm'
    ds.depth_aqc.attrs['units'] = 'm'
    ds.r.attrs['long_name'] = 'realised growth rate'
    ds.r.attrs['units'] = 'days^-1'
    ds.stdev_r.attrs['long_name'] = 'standard deviation in realised growth rate'
    ds.stdev_r.attrs['units'] = 'days^-1'
    ds.stdev_r_norm.attrs['long_name'] = 'normalised standard deviation in realised growth rate'
    ds.stdev_r_norm.attrs['units'] = '-'
    ds.stdev_chl.attrs['long_name'] = 'standard deviation in chlorophyll concentration'
    ds.stdev_chl.attrs['units'] = 'mg m^-3'
    ds.stdev_chl_norm.attrs['long_name'] = 'coefficient of variation in chlorophyll concentration'
    ds.stdev_chl_norm.attrs['units'] = '-'
    # add global attributes
    ds.attrs["search radius_r"] = 'search radius for computing standard deviations = '+str(rad)+' m'
    # write the output
    ds.to_netcdf(fname_out)
    
if __name__ == "__main__":
    
    react_params=params.reactions01 # select which reactions params to use
    physics_params=params.physics02 # select which physics params to use
    
    ps = [1.e-3, 1.e-7]
    for p in ps:
        do_the_diags(react_params,physics_params,p)
    
    
