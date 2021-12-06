import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
import math
ion()
from plot_croco_lib import *

Qswmaxs=[250,800]
tau_means=[0,0.05]

for tau_mean in tau_means:
    for Qswmax in Qswmaxs:
        
        fig, ax = plt.subplots(figsize=(10,20))

        ax = [subplot(3,1,i+1) for i in range(3)]
        
        for n in range(3):
            ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
        
        amplitudes = [0.01, 0.02,0.03, 0.04]
        for amplitude in amplitudes:
            crocofile = 'mean'+str(tau_mean)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
            
            time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
            z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
            # get kappa onto the rho axis so each value will have a cell depth
            kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
            kappa_croco_surf=get_Cs_eulerian(time_croco,z,zw,kappa_croco_r,z_therm_croco)
            
            ax[1].plot(time_croco,kappa_croco_surf,label=str(amplitude))
            ax[2].plot(time_croco,z_therm_croco,label=str(amplitude))
        
        for n in range(3):
            ax[n].set_xlim(1,14)
            ax[n].set_xticks(range(0,15))
            
        # sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
        ax[0].plot(time_croco,s_flux,'k')
        ax[0].set_xticklabels([])
        ax[0].set_ylabel('Surface radiation (W m$^{-2}$)', fontsize=15,)
        ax[0].set_ylim(0,800)
        
        ax[1].set_ylabel("$\kappa_s$ (m$^2$ s$^{-1}$)", fontsize=15)
        ax[1].set_xticklabels([])
        ax[1].set_ylim(0, 0.003)
        
        ax[2].set_ylabel("$\ell$ (m)", fontsize=15)
        ax[2].set_xlabel("Time (days)", fontsize=15)
        ax[2].set_ylim(10, 25)
        
        ax[2].legend(fontsize=10, loc="upper left",title="$\\tau^{ac0}$ (N m$^{-2}$)")
            
        plt.savefig('plot_kappa_mean'+str(tau_mean)+'_flx'+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
