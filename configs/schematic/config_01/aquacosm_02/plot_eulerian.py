import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
ion()

def plot_C(ax,n,time,z,carbon,label,max_chl):
    #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
    img = ax[n].scatter(time, z, c=carbon,
                        cmap=cm.viridis, linewidth=0, s=40)
    img.set_clim(0, max_chl)
    ax[n].set_ylim(50, 0)
    ax[n].set_xlim(0,21)
    ax[n].set_xticks(range(0,22))
    ax[n].set_ylabel('Depth (m)', fontsize=15,)
    if n==0:
        # ax[n].set_xlabel('Time (days)', fontsize=15,)
        cbarax = gcf().add_axes((0.82, 0.55, 0.015, 0.5)) # [left, bottom, width, height] 
        cbar = colorbar(img, cbarax)
        cbar.set_label("Chlorophyll  (mg m$^{-3}$)", fontsize=15)
    # else:
    #     ax[n].set_xticklabels([])
    
    ax[n].set_xticklabels([])


figure(figsize=(10,5))
#ax = [subplot(4,1,i+1) for i in range(4)]
ax = [subplot(1,1,i+1) for i in range(1)]   

# max values to plot
max_chl=25

# plot the eulerian data
eulfile='eulerian_r1.0_c0.0_a0.5_l5.0_mld50_kappa0.0001_dt5.0.nc'
time_eul,z_eul,chl_eul,chl_eul_avg=get_eul_output(eulfile)

Nt_eul,Nz_eul=shape(chl_eul)
# repeat time along the z dimension for plotting
time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
# repeat depth along the time dimension for plotting
z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)

plot_C(ax,0,time_eul,z_eul,chl_eul,'Eulerian',max_chl)
    
# plot the average chlorophyll
ts=gcf().add_axes((0., 0.55-0.55*1, 0.8, 0.5))
ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=4, label='Eulerian')
ts.set_ylabel('average Chl (mg m$^{-3}$)', fontsize=15)
ts.set_xlabel('Time (days)', fontsize=15)
ts.set_ylim(0, max_chl)
ts.set_xlim(0,21)
ts.set_xticks(range(0,22))
ts.grid(linestyle=':', linewidth=0.5)

plt.savefig(Path(eulfile).stem+'.jpg',dpi=500,bbox_inches = 'tight')
