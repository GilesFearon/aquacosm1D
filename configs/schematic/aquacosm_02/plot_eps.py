import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
import math
ion()


def get_r(time,Cs):
    r = np.zeros(len(time))
    r[:] = np.NaN
    
    dCs_dt = (Cs[1:]-Cs[0:-1])/(time[1:]-time[0:-1])
    
    r[0:-1] = dCs_dt / Cs[1:]
    
    return r

def get_r_reactions(time,z,C,reaction):
    r = np.zeros(len(time))
    r[:] = np.NaN
    
    crocodir='../physics/'
    crocofilename="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
    crocofile=crocodir+crocofilename                
    dt = 5             
    wc = water_column_netcdf(DatasetName=crocofile, max_depth=mld)
    React = set_up_reaction(wc, dt, reaction, 
                               LightDecay = 5.,
                               BasePhotoRate = 1.,
                               RespirationRate = 0.1)
    Nt,Npts=shape(C)
    
    Nscalars=1

    Tracers      = Aquacosm1D_Particles(
            zeros((Npts, Nscalars+2), dtype='float64')
            )
    Tracers[:,0] = np.arange(Npts)
    Tracers[:,1] = z
    #Tracers[:,2] += 1
    
    for ii in range(len(time)):
        React.wc.set_current_time(time[ii])
        Tracers[:,2]=C[ii,:]
        RRates=React.current_model(Tracers, React.wc, time[ii])
        r[ii]=np.abs(np.mean(RRates[:,0]/Tracers[:,2])) # abs to handle when reactions go negative
        
    return r
    
def get_eul_output(eulfile):
    data_eul=Dataset(eulfile)
    # Chl_C=np.float(data_eul.Chl_C) # ratio of chlorophyll to carbon [mgChl/mgC]
    chl_eul=data_eul.variables['chl'][:,1:-1] 
    z_eul=data_eul.variables['depth'][1:-1] # cropping end points used in eulerian code as these are only used as boundary conditions 
    time_eul=data_eul.variables['time'][:]/3600/24 # time in days
    time_eul=time_eul-time_eul[0]
    chl_eul_avg = np.mean(chl_eul, axis=1)
    Nt_eul,Nz_eul=shape(chl_eul)
    # # repeat time along the z dimension for plotting
    # time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # # repeat depth along the time dimension for plotting
    # z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    data_eul.close()
    return time_eul,z_eul,chl_eul,chl_eul_avg

fig, ax = plt.subplots(figsize=(10,5))

ax = [subplot(2,1,i+1) for i in range(2)]

# physical inputs to loop through for sensitivity tests
mlds = [20, 50]
kappas = [0.0001, 0.001] #[0.0001,0.001,0.01]  
max_photo = '1.0'
reaction = Sverdrup_incl_K #'Sverdrup_incl_K'

for n,mld in enumerate(mlds):
    ax[n].set_position(  (0., 0.55-0.6*n, 0.8, 0.5))
    for kappa in kappas: 
        
        #eulfile='eulerian_'+reaction+'_r'+max_photo+'_mld'+str(mld)+"_kappa"+str(kappa)+".nc"
        eulfile='eulerian_'+reaction.__name__+'_r'+max_photo+'_mld'+str(mld)+"_kappa"+str(kappa)+".nc"
        time_eul,z_eul,chl_eul,chl_eul_avg=get_eul_output(eulfile)
        
        # get r from the model output
        # not bad for the 20 m case, but it is not constant when it should be
        # r=get_r(time_eul,chl_eul_avg)/86400 # in s^-1
        
        # get r directly from the reactions function        
        r = get_r_reactions(time_eul*3600*24,z_eul,chl_eul,reaction)
                
        # compute epsilon
        ell = mld
        # ell = 20 # should we not compute based on a euphtic zone? This will also influence r, which should only get meaned over this depth
        
        eps=r*np.square(ell)/kappa
        
        ax[n].plot(time_eul,eps,label='$\kappa = $'+str(kappa)+' $m^2 s^{-1}$')

    ax[n].set_xlim(0,21)
    
    ax[n].set_ylabel("$\epsilon$", fontsize=15)
    ax[n].set_yscale('log') 
    
    if n==1:
        ax[n].set_xlabel("Time (days)", fontsize=15)
    else:
        ax[n].set_xticklabels([])
    
    ax[n].set_title('$\ell = $'+str(mld)+ 'm')
    
    ax[n].set_ylim(1.e-4, 1.e2)
        
    ax[n].legend(fontsize=10, loc="upper left")
    
plt.savefig('plot_eps_'+reaction.__name__+'_r'+max_photo+'.jpg',dpi=500,bbox_inches = 'tight')
