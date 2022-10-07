import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import xarray as xr
from pathlib import Path
import params # params.py must be in this directory

react_params=params.reactions01 # select which reactions params to use
physics_params=params.physics02 # select which physics params to use

# seed(1234567) moving below to reset seed for each run

#------------------------------------------------------------
dt        = 5. # time step in seconds
Ndays     = 21 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(0.5*3600 / dt) #store the particles every Nshow time steps
Nconsole  = int(6*3600 / dt) # frequency of writing to the console
Npts      = 200  #number of particles
Nscalars  = 1    #number of scalars carried by each particle

# aquacosm settings to loop through for sensitivity tests
ps = [1.e-3, 1.e-7]
 
for p in ps:

    seed(1234567)
    
    crocodir='../physics/'
    crocofilename="mean"+str(physics_params.mean_tau)+"_mld"+str(physics_params.mld)+"_amp"+str(physics_params.amplitude)+"_flx"+str(physics_params.Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename                
               
    wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
    
    Diffuse = set_up_diffusion(Npts, Nscalars,
                               radius=10.,
                               p=p, #1.e-3, ###################
                               dt=dt,
                               wc=wc)
    
    Transport = set_up_transport(wc, dt, SODE='Milstein')
    
    # Reactions
    # BioShading_onlyC
    React = set_up_reaction(wc, dt, BioShading_onlyC,
                        LightDecay=react_params.LightDecay,
                        AlphaEpsilon=react_params.AlphaEpsilon,
                        MaxPhotoRate = react_params.MaxPhotoRate, 
                        BasalMetabolism = react_params.BasalMetabolism,
                        Chl_C = react_params.Chl_C,
                        CrowdingMortality = react_params.CrowdingMortality,
                        CrowdingHalfSaturation = react_params.CrowdingHalfSaturation,
                        Chl_light_abs = react_params.Chl_light_abs)
            
    Particles = create_particles(Npts, Nscalars, wc)
    # Here's where we initialise the chlorophyll value for the particles
    data_croco=Dataset(crocofile)
    tpas=data_croco.variables['tpas'][0,:,0,0]
    temp=data_croco.variables['temp'][0,:,0,0]
    zt=data_croco.variables['deptht'][:]
    data_croco.close()
    # create a constant chlorophyll ini over surface layer
    temp_thermocline=11
    z_therm=interp1d(temp,zt,kind='linear')(temp_thermocline)
    chl_ini=np.zeros(np.shape(zt))+1e-3 # mg/m3
    chl_ini[zt<z_therm]=1
    # extend the zt and chl_ini arrays for interpolation to particle locations near boundaries
    chl_ini=concatenate(([chl_ini[0]], chl_ini, [chl_ini[-1]]))
    zt=concatenate(([-10], zt,[9999]))
    chl_ini = interp1d(zt,chl_ini,kind='linear')(Particles[:,1])
    # convert Chl to C using fixed ratio
    c_ini=chl_ini/React.Chl_C  
    Particles[:, 2] = c_ini
    
    fname_out='aquacosm_p'+"{0:1.0e}".format(Diffuse.p)+'_'+react_params.Name+'_'+physics_params.Name+'.nc'
    print('working on: ', fname_out)
    
    Pstore    = []   #list that stores the particles every Nstore time steps
    Tstore    = []   #output time stamps
    starttime = wc.times[0]
    wc.set_current_time(starttime)
    for loop in range(Nloops):
        #-------------------------
        if loop % Nstore == 0:
            sort_by_depth(Particles)
            Pstore.append(Particles.copy())
            Tstore.append(wc.get_current_time())
            
        if loop % Nconsole == 0:
            print('Days elapsed: ', loop*dt/(24*3600.), 'Mean C',
                  mean(Particles[:,2]))
        
        Diffuse(Particles)
        React(Particles)
        Transport(Particles)
    
        wc.increment_current_time(dt)
        
        # here is where we repeat the input data if we reach the end of the input time
        # I'm not considering this for now
        #
        # try:
        #     wc.increment_current_time(dt)
        # except ValueError: 
        #     #Save partial results
        #     Pst = array(Pstore)
        #     Pst.tofile("aquacosm_analytical.pydat".format(Diffuse.p))
        #     wc.set_current_time(starttime)
    
    Pstore = array(Pstore)
    
    # write a netcdf output file
    ds = xr.Dataset(data_vars=dict(depth=(["time","depth_rank"],Pstore[:,:,1]),
                                           chl=(["time","depth_rank"],Pstore[:,:,2]*React.Chl_C), # converting to Chlorophyll concentration for output
                                           aqc_id=(["time","depth_rank"],Pstore[:,:,0]),
                                           ),
                                           coords=dict(depth_rank=(["depth_rank"],Pstore[0,:,0]),time=(["time"],Tstore)))
    ds.time.attrs['units'] = 'seconds since start of simulation' # not a real time period as an analytical experiment
    ds.depth_rank.attrs['long_name'] = 'depth rank of aquacosms'
    ds.depth_rank.attrs['units'] = '-'
    ds.aqc_id.attrs['long_name'] = 'aquacosm unique identifier'
    ds.aqc_id.attrs['units'] = '-'
    ds.chl.attrs['long_name'] = 'chlorophyll concentration'
    ds.chl.attrs['units'] = 'mg m^-3'
    ds.depth.attrs['long_name'] = 'depth of aquacosm'
    ds.depth.attrs['units'] = 'm'
    # add global attributes telling us what the input parameters were
    ds.attrs["couping_parameter_p"] = "{:.1e}".format(Diffuse.p)
    ds.attrs["Chl_C"] = str(React.Chl_C)
    # write the output
    ds.to_netcdf(fname_out)
