import sys
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_croco_lib import *
ion()

def plot_C(ax,n,time,z,data,label,cmap,min_value,max_value):
    #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
    img = ax[n].scatter(time, z, c=data,
                        cmap=cmap, linewidth=0, s=40)
    img.set_clim(min_value, max_value)
    ax[n].set_ylim(50, 0)
    ax[n].set_xlim(0,14)
    ax[n].set_xticks(range(0,15))
    ax[n].set_ylabel('Depth (m)', fontsize=15,)
    if n==2:
        ax[n].set_xlabel('Time (days)', fontsize=15,)
    else:
        ax[n].set_xticklabels([])
    
    cbarax = gcf().add_axes((0.82, 0.55-0.55*n, 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    # cbar.set_label(label, fontsize=15)
    
    ax[n].text(0.2, 45, label,
               fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    
def do_the_plot(mld,amplitude,Qswmax):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(3,1,i+1) for i in range(3)]   
    
    # get the croco input data
    crocodir='./'
    crocofilename="mean0_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename   
    time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
    z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
    Nt_croco,Nz_croco=shape(temp_croco)
    # repeat time along the z dimension for plotting
    time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
    # repeat depth along the time dimension for plotting
    z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0)
    
    plot_C(ax,0,time_croco,z_croco,u_croco,'u (m s$^{-1}$)',cm.RdBu,-0.4,0.4)
    plot_C(ax,1,time_croco,z_croco,temp_croco,'temperature ('+u"\N{DEGREE SIGN}"+' C)',cm.jet,10,18)
    plot_C(ax,2,time_croco,z_croco,tpas_croco,'tracer (-)',cm.viridis,0,1)
    
    # add the thermocline
    for n in range(0,3):
        cnt=ax[n].contour(time_croco,z_croco,temp_croco,[11],colors='w',linewidths=2.5,linestyles='dashed')
        
    # sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
    # sx.plot(time_croco[:,0],s_flux)
    # sx.set_xticklabels([])
    # sx.set_ylabel('Surface heat flux (W m$^{-2}$)', fontsize=15,)
    # sx.set_xlim(0,21)
    # sx.set_xticks(range(0,22))
    
    plt.savefig('plot_croco_output_'+Path(crocofile).stem+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    amplitudes = [0.03] #[0, 0.01, 0.02, 0.03, 0.04]
    mlds = [10] #[10, 25]   
    Qswmaxs = [800]
        
    for amplitude in amplitudes:
        for mld in mlds:
            for Qswmax in Qswmaxs:
                do_the_plot(mld,amplitude,Qswmax)
    
    
