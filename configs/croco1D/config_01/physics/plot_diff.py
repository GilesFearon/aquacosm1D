import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
import math
ion()

def get_kappa_surface(time,z,zw,temp,kappa,temp_thermocline):
    kappa_s=np.zeros(len(time))
    z_therm=np.zeros(len(time))   
    
    z_thickness=zw[1:]-zw[0:-1];
    
    for t in range(0,len(time)):
        #could put a check here to make sure temp_thermocline is in the range of temp[t,:]
        z_therm[t]=interp1d(temp[t,:],z,kind='linear')(temp_thermocline)
        # find the z indices above and below the interpolated depth
        indx_above=np.where(zw<z_therm[t])[-1][-1]
        indx_below=indx_above+1
        # compute the weights to apply to the indices either side of the depth of max strat
        weight_above=(z[indx_above]-z_therm[t])/(z[indx_above]-z[indx_below]);
        weight_below=1-weight_above;
        #
        # integrate kappa from the surface to z_therm
        for k in range(0,indx_above): # using range means we don't include indx_above in the loop, which is what we want
            kappa_s[t]=kappa_s[t] + kappa[t,k]*z_thickness[k]
        # add half of the layer above
        kappa_s[t]=kappa_s[t]+0.5*kappa[t,indx_above]*z_thickness[indx_above]
        # add the mean of the layer above and below, weighted by where the
        # thermocline is relative to the two layers
        kappa_s[t]=kappa_s[t]+weight_above*np.mean(kappa[t,indx_above:indx_below])*np.mean(z_thickness[indx_above:indx_below])
    
    kappa_s[kappa_s == 0] = np.NaN 
    z_therm[z_therm == 0] = np.NaN   

    kappa_s = kappa_s/z_therm # now the mean over the surface layer    
         
    return kappa_s, z_therm   
    
if __name__ == "__main__":
    
    fig, ax = plt.subplots(figsize=(10,15))
    
    ax = [subplot(3,1,i+1) for i in range(3)]
    
    for n in range(3):
        ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    
    amplitudes = [0.01, 0.02,0.03, 0.04]
    for amplitude in amplitudes:
        crocofile = 'mean0_mld10_amp'+str(amplitude)+'_flx800_lat30_T016_hmax50.nc'
        
        # get the croco output
        data_croco=Dataset(crocofile)
        # tpas_croco=data_croco.variables['tpas'][:,:,0,0]
        temp_croco=data_croco.variables['temp'][:,:,0,0]
        kappa_croco=data_croco.variables['difvho'][:,:,0,0]
        z=data_croco.variables['deptht'][:]
        zw=data_croco.variables['depthw'][:]
        time_croco=data_croco.variables['time_counter'][:]/3600/24 # time in days
        kappa_croco_surf,z_therm_croco=get_kappa_surface(time_croco,z,zw,temp_croco,kappa_croco,11)
        # compute epsilon
        r=0.5 # growth rate days^-1
        eps=r/86400*np.square(z_therm_croco)/kappa_croco_surf
        data_croco.close()
        
        h=50/200
        dt=5 #s
        ps = [1e-3, 1e-5, 1e-7]
        for i, p in enumerate(ps):
            q = p/(4*np.pi*kappa_croco_surf*dt)*exp(-1*math.pow(h,2)/(4*kappa_croco_surf*dt))
            D = math.pow(h,2)/dt*q
            ax[i].plot(time_croco,D,label=str(amplitude))
            ax[i].plot([0,8],[1e-9,1e-9],color='k',linestyle='dotted', linewidth=0.5)
            if i<2:
                ax[i].set_xticklabels([])
            else:
                ax[i].set_xlabel("Time (days)", fontsize=15)
                ax[i].legend(fontsize=10, loc="upper left",title="$tau^{ac0}$ (N m$^{-2}$)")
            ax[i].set_xlim(1,8)
            ax[i].set_yscale('log')
            ax[i].set_ylim(1e-18, 1e-3)
            ax[i].set_ylabel("D (m$^2$ s$^{-1}$)\np = "+"{:.1e}".format(p), fontsize=15)
        
        
    plt.savefig('plot_diff.jpg',dpi=500,bbox_inches = 'tight')
    
