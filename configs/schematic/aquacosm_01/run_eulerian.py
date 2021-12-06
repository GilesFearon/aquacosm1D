import sys
sys.path.insert(0, '../../../aquacosm1D_lib')
from aquacosm1D_reactions   import set_up_reaction, NoReactions, BioShading_onlyC, Sverdrup, Sverdrup_incl_K
from aquacosm1D_watercolumn import water_column, water_column_netcdf
from eulerian1D_diffusion import set_up_eulerian_diffusion
from aquacosm1D_utilities import Aquacosm1D_Particles
from pylab import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.interpolate import splev
import xarray as xr

ion()

seed(1234567)

#------------------------------------------------------------
dt        = 5. # time step in seconds
Ndays     = 7 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(0.5*3600 / dt) #store the particles every Nshow time steps
Nconsole  = int(6*3600 / dt) # frequency of writing to the console
Nscalars  = 1    #number of scalars carried by each particle

# physical inputs to loop through for sensitivity tests
mlds = [50] #[20,50]
kappas = [0.0001] #[0.0001,0.001,0.01]  

for mld in mlds:
    for kappa in kappas:          
        crocodir='../physics/'
        crocofilename="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
        crocofile=crocodir+crocofilename                
                        
        wc = water_column_netcdf(DatasetName=crocofile, max_depth=mld)
        
        Diffuse = set_up_eulerian_diffusion(wc,dt,Nscalars)
                        
        # React = set_up_reaction(wc, dt, BioShading_onlyC,
        #                         LightDecay=10.,
        #                         MaxPhotoRate = 2., 
        #                         BasalMetabolism = 0.16,
        #                         Chl_C = 0.017,
        #                         CrowdingMortality = 0.)
        #
        # React = set_up_reaction(wc, dt, NoReactions)
        # React.Chl_C = 0.017 # just adding this here to be consistent with what we write out with BioShading_onlyC
        #
        # React = set_up_reaction(wc, dt, Sverdrup, 
        #                         LightDecay = 5.,
        #                         BasePhotoRate = 1.,
        #                         RespirationRate = 0.1)
        # React.Chl_C = 1. # Not applicable. Just adding this here for compatibility with what we write out with BioShading_onlyC
        
        React = set_up_reaction(wc, dt, Sverdrup_incl_K, 
                                    LightDecay = 5.,
                                    BasePhotoRate = 10.,
                                    RespirationRate = 0.1,
                                    CarryingCapacity = 20)
        React.Chl_C = 1. # Not applicable. Just adding this here for compatibility with what we write out with BioShading_onlyC
            
        # Here's where we initialise the chlorophyll       
        # create the eulerian tracer array in the same format as the aquacosms so we can 
        # use the same reactions library as used in the aquacosms
        Npts=len(Diffuse.zc) 
        Tracers      = Aquacosm1D_Particles(
                zeros((Npts, Nscalars+2), dtype='float64')
                )
        Tracers[:,0] = np.arange(Npts)
        Tracers[:,1] = Diffuse.zc
        Tracers[:,2] += 5
        
        fname_out='eulerian_'+type(React.current_model).__name__+'_r'+str(React.BasePhotoRate*(60.*60.*24.))+"_mld"+str(mld)+"_kappa"+str(kappa)+".nc"
        print('working on: ', fname_out)
        
        Cstore = []
        Tstore    = []   #output time stamps
        starttime = wc.times[0]
        wc.set_current_time(starttime)
        for loop in range(Nloops):
            #-------------------------
            
            if loop%Nstore == 0:
                Cstore.append(Tracers.copy())
                Tstore.append(wc.get_current_time())
                
            if loop % Nconsole == 0:        
                print("Days elapsed:", loop*dt/(24*3600.),
                      'Mean C:', mean(Tracers[:,2]))
                
            Diffuse(Tracers)
            React(Tracers)
            
            wc.increment_current_time(dt)
        
        Cstore = array(Cstore)
        Tstore = array(Tstore)
        
        # write a netcdf output file
        ds = xr.Dataset(data_vars=dict(chl=(["time","depth"],Cstore[:,:,2]*React.Chl_C)),
                                               coords=dict(depth=(["depth"],Diffuse.zc),time=(["time"],Tstore)))
        ds.time.attrs['units'] = 'seconds since start of simulation' # not a real time period as an analytical experiment
        ds.chl.attrs['long_name'] = 'chlorophyll concentration'
        ds.chl.attrs['units'] = 'mg m^-3'
        ds.depth.attrs['long_name'] = 'depth positive downward'
        ds.depth.attrs['units'] = 'm'
        # add global attributes telling us what the input parameters were
        ds.attrs["Chl_C"] = str(React.Chl_C)
        # write the output
        ds.to_netcdf(fname_out)
        
