import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import xarray as xr
from pathlib import Path

# seed(1234567) moving below to reset seed for each run

#------------------------------------------------------------
dt        = 5. # time step in seconds
Ndays     = 21 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(0.5*3600 / dt) #store the particles every Nshow time steps
Nconsole  = int(6*3600 / dt) # frequency of writing to the console
# Npts      = 200  #number of particles
Nscalars  = 1    #number of scalars carried by each particle

# physical inputs to loop through for sensitivity tests
# (corresponding to CROCO 1D runs)
mlds = [20, 50] #[20,50]
kappas = [0.001, 0.0001] #[0.0001,0.001,0.01] 

# aquacosm settings to loop through for sensitivity tests
ps = [1e-3,1e-7]

for mld in mlds:
    Npts = mld*4  #scale number of particles with depth
    for kappa in kappas:
        for p in ps:
        
            seed(1234567)
            
            crocodir='../physics/'
            crocofilename="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
            crocofile=crocodir+crocofilename                
                       
            wc = water_column_netcdf(DatasetName=crocofile, max_depth=mld)
            
            Diffuse = set_up_diffusion(Npts, Nscalars,
                                       radius=10.,
                                       p=p, #1.e-3, ###################
                                       dt=dt,
                                       wc=wc)
            
            Transport = set_up_transport(wc, dt, SODE='Milstein')
            
            # Reactions
            # BioShading_onlyC
            React = set_up_reaction(wc, dt, BioShading_onlyC,
                                LightDecay=5.,
                                MaxPhotoRate = 1., 
                                BasalMetabolism = 0.16,
                                Chl_C = 0.017,
                                CrowdingMortality = 0.5,
                                Chl_light_abs = 0.)
                    
            Particles = create_particles(Npts, Nscalars, wc)
            # Here's where we initialise the chlorophyll value for the particles
            chl_ini = Particles[:, 2] + 1 # mg/m3 (Particles[:, 2] is initialised as zeros)
            # convert Chl to C using fixed ratio
            c_ini=chl_ini/React.Chl_C  
            Particles[:, 2] = c_ini
            
            fname_out='aquacosm_p'+"{0:1.0e}".format(Diffuse.p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mld'+str(mld)+'_kappa'+str(kappa)+'_dt'+str(dt)+'.nc'
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
            #Pstore.tofile("quacosm_analytical.pydat".format(Diffuse.p))
            
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
