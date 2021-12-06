import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
ion()

def plot_C(ax,n,time,z,carbon,label,t,max_chl):
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
    
    # ax[n].text(0.2, 45, label,
    #            fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    # ax[n].plot([t, t], [0, 50], '--k',
    #                     linewidth=1)

def do_the_plot(mld,kappa,reaction):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(1,1,i+1) for i in range(1)]   
    
    # time snapshot to plot
    t=5.75 #days
    
    # max values to plot
    if reaction=='NoReactions':
        max_chl = 1
    elif reaction=='Sverdrup_incl_K':
        max_chl = 20
    else:
        max_chl = 10
        
    # plot the eulerian data
    eulfile='eulerian_'+reaction+'_'+'mean0_amp'+str(amplitude)+'_mld'+str(mld)+'_flx250.nc'
    time_eul,z_eul,chl_eul=get_eul_output(eulfile)
    Nt_eul,Nz_eul=shape(chl_eul)
    # repeat time along the z dimension for plotting
    time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # repeat depth along the time dimension for plotting
    z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    tindx_eul = (np.abs(time_eul[:,0] - t)).argmin()
    plot_C(ax,0,time_eul,z_eul,chl_eul,'Eulerian',t,max_chl)
   
    plt.savefig(Path(eulfile).stem+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":

    amplitudes = [0.01] #[0, 0.01, 0.02, 0.03, 0.04]
    mlds = [10] #[10, 25]   
    reactions = ['Sverdrup_incl_K'] #['Sverdrup','Sverdrup_incl_K']
        
    for amplitude in amplitudes:
        for mld in mlds:
            for reaction in reactions:
                do_the_plot(mld,amplitude,reaction)
    
    
