import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
ion()

Qswmaxs=[250]
tau_means=[0]
reactions = [Sverdrup_incl_K]

for tau_mean in tau_means:
    for Qswmax in Qswmaxs:
        for reaction in reactions:
        
            fig, ax = plt.subplots(figsize=(10,20))
    
            ax = [subplot(3,1,i+1) for i in range(3)]
            
            for n in range(3):
                ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
            
            amplitudes = [0.01, 0.03]
            for amplitude in amplitudes:
                crocodir='../physics/'
                crocofilename = 'mean'+str(tau_mean)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
                crocofile=crocodir+crocofilename 
                
                time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
                z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
                # get kappa onto the rho axis so each value will have a cell depth
                kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
                
                # get the eulerian data
                eulfile='eulerian_'+reaction.__name__+'_mean0_amp'+str(amplitude)+'_mld10_flx'+str(Qswmax)+'.nc'
                time_eul,z_eul,chl_eul=get_eul_output(eulfile)
                # interpolate z_therm onto eulerian time axis (just in case different)
                z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
                
                # interpolate kappa_croco_r onto eulerian time axis
                kappa_eul_r=np.zeros_like(chl_eul)
                for ii,zi in enumerate(z):
                    kappa_eul_r[:,ii]=interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([kappa_croco_r[0,ii]],kappa_croco_r[:,ii])),kind='linear')([time_eul])
                
                # get r directly from the reactions function        
                r,z_euphotic = get_r_reactions(crocofile,time_eul*3600*24,z,zw,chl_eul,z_therm_eul,reaction)
                
                kappa_eul_euphotic = get_Cs_eulerian(time_eul,z,zw,kappa_eul_r,z_euphotic)
                
                # compute epsilon
                eps=r*np.square(z_euphotic)/kappa_eul_euphotic
                
                ax[0].plot(time_eul,kappa_eul_euphotic,label=str(amplitude))
                ax[1].plot(time_eul,z_euphotic,label=str(amplitude))
                ax[2].plot(time_eul,eps,label=str(amplitude))
            
            for n in range(3):
                ax[n].set_xlim(1,21)
                ax[n].set_xticks(range(0,22))
            
            
            # sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
            # sx.plot(time_croco,s_flux,'k')
            # sx.set_xticklabels([])
            # sx.set_ylabel('Surface radiation (W m$^{-2}$)', fontsize=15,)
            # sx.set_ylim(0,800)
            
            ax[0].set_ylabel("$\kappa_s$ (m$^2$ s$^{-1}$)", fontsize=15)
            ax[0].set_xticklabels([])
            ax[0].set_ylim(0, 0.001)
            
            ax[1].set_ylabel("$\ell$ (m)", fontsize=15)
            ax[1].set_xticklabels([])
            ax[1].set_ylim(10, 25)
            
            ax[2].set_ylabel("$\epsilon$ (-)", fontsize=15)
            ax[2].set_xlabel("Time (days)", fontsize=15)
            ax[2].set_ylim(0, 10)
            
            ax[0].legend(fontsize=10, loc="upper left",title="$\\tau^{ac0}$ (N m$^{-2}$)")
                
            plt.savefig('plot_timeseries_'+reaction.__name__+'_mean'+str(tau_mean)+'_mld10_flx'+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
