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

# def plot_C(ax,n,time,z,carbon,label):
#     #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
#     # ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
#     img = ax[n].scatter(time, z, c=carbon,
#                         cmap=cm.viridis, linewidth=0, s=40)
#     img.set_clim(0, 1)
#     ax[n].set_ylim(50, 0)
#     # ax[n].set_xlim(0,14)
#     # ax[n].set_xticks(range(0,15))
#     ax[n].set_ylabel('Depth (m)', fontsize=15,)
#     if n==2:
#         # ax[n].set_xlabel('Time (days)', fontsize=15,)
#         cbarax = gcf().add_axes((0.82, -0.525, 0.015, 1.)) # [left, bottom, width, height] 
#         cbar = colorbar(img, cbarax)
#         cbar.set_label("tracer (-)", fontsize=15)
#     # else:
#     #     ax[n].set_xticklabels([])
    
#     #ax[n].set_xticklabels([])
    
#     ax[n].text(0.2, 45, label,
#                fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

def plot_C(ax,n,time,z,data,label,cmap,min_value,max_value):
    #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
    img = ax[n].scatter(time, z, c=data,
                        cmap=cmap, linewidth=0, s=40)
    img.set_clim(min_value, max_value)
    ax[n].set_ylim(50, 0)
    # ax[n].set_xlim(0,14)
    # ax[n].set_xticks(range(0,15))
    ax[n].set_ylabel('Depth (m)', fontsize=15,)
    # if n==2:
    #     ax[n].set_xlabel('Time (days)', fontsize=15,)
    # else:
    #     ax[n].set_xticklabels([])
    
    cbarax = gcf().add_axes((0.82, 0.55-0.55*n, 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    # cbar.set_label(label, fontsize=15)
    
    ax[n].text(0.2, 45, label,
               fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    
Qswmaxs=[800]
tau_means=[0]

for tau_mean in tau_means:
    for Qswmax in Qswmaxs:
        
        figure(figsize=(10,5))
        
        ax = [subplot(4,1,i+1) for i in range(4)]
        
        for n in range(4):
            ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
            # ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
        
        # amplitudes = [0.01, 0.03]
        # for ii,amplitude in enumerate(amplitudes):
        amplitude = 0.03
            
        crocofile = 'mean'+str(tau_mean)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
        time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
        z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
        # get kappa onto the rho axis so each value will have a cell depth
        kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
        kappa_croco_surf=get_Cs_eulerian(time_croco,z,zw,kappa_croco_r,z_therm_croco)
        
        Nt_croco,Nz_croco=shape(temp_croco)
        # repeat time along the z dimension for plotting
        time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
        # repeat depth along the time dimension for plotting
        z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0)
        
        # plot the passive tracer
        plot_C(ax,1,time_croco,z_croco,temp_croco,'temperature ('+u"\N{DEGREE SIGN}"+' C)',cm.jet,10,18)
        plot_C(ax,2,time_croco,z_croco,tpas_croco,'tracer (-)',cm.viridis,0,1)
         # add the thermocline
        for n in range(1,3):
            ax[n].contour(time_croco,z_croco,temp_croco,[11],colors='w',linewidths=2.5,linestyles='dashed')
        
        # plot the surface averaged kappa
        ax[3].plot(time_croco[:,0],kappa_croco_surf,label=str(amplitude))
        
        for n in range(4):
            ax[n].set_xlim(1,14)
            ax[n].set_xticks(range(0,15))
            if n==3:
                ax[n].set_xlabel("Time (days)", fontsize=15)
            else:
                ax[n].set_xticklabels([])
        
        ax[0].plot(time_croco[:,0],s_flux,'k')
        ax[0].set_ylabel('$Q_s$ (W m$^{-2}$)', fontsize=15,)
        ax[0].set_ylim(0,800)
        
        ax[3].set_ylabel("$\kappa_s$ (m$^2$ s$^{-1}$)", fontsize=15)
        ax[3].set_ylim(0, 0.0025)
        # ax[3].legend(fontsize=10, loc="upper left",title="$\\tau^{ac0}$ (N m$^{-2}$)")
            
        plt.savefig('plot_tpas_temp_kappa_mean'+str(tau_mean)+'_flx'+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
