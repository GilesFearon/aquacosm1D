import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
ion()

def get_eul_output(eulfile):
    data_eul=Dataset(eulfile)
    # Chl_C=np.float(data_eul.Chl_C) # ratio of chlorophyll to carbon [mgChl/mgC]
    chl_eul=data_eul.variables['chl'][:,1:-1] 
    z_eul=data_eul.variables['depth'][1:-1] # cropping end points used in eulerian code as boundary conditions - produces same grid as croco by definition 
    time_eul=data_eul.variables['time'][:]/3600/24 # time in days
    time_eul=time_eul-time_eul[0]
    chl_eul_avg = np.mean(chl_eul, axis=1)
    Nt_eul,Nz_eul=shape(chl_eul)
    # repeat time along the z dimension for plotting
    time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # repeat depth along the time dimension for plotting
    z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    data_eul.close()
    return time_eul,z_eul,chl_eul,chl_eul_avg

def plot_C(ax,n,time,z,carbon,label,t,max_chl):
    #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
    img = ax[n].scatter(time, z, c=carbon,
                        cmap=cm.viridis, linewidth=0, s=40)
    img.set_clim(0, max_chl)
    ax[n].set_ylim(50, 0)
    ax[n].set_xlim(0,7)
    ax[n].set_xticks(range(0,7))
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

def do_the_plot(mld,kappa,reaction,max_photo):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(1,1,i+1) for i in range(1)]   
    
    # time snapshot to plot
    t=5.75 #days
    
    # max values to plot
    if reaction=='Sverdrup_incl_K':
        max_chl = 20
    else:
        max_chl = 30
        
    # plot the eulerian data
    eulfile='eulerian_'+reaction+'_r'+max_photo+'_mld'+str(mld)+"_kappa"+str(kappa)+".nc"
    time_eul,z_eul,chl_eul,chl_eul_avg=get_eul_output(eulfile)
    tindx_eul = (np.abs(time_eul[:,0] - t)).argmin()
    plot_C(ax,0,time_eul,z_eul,chl_eul,'Eulerian',t,max_chl)
        
    # compare the chlorophyll averaged over the surface layer
    # ts=gcf().add_axes((0.4, -0.8, 0.45, 0.6))
    ts=gcf().add_axes((0., 0.55-0.55*1, 0.8, 0.5))
    ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=4, label='Eulerian')
    ts.set_ylabel('average surafce Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    ts.set_ylim(4, max_chl)
    ts.set_xlim(0,7)
    ts.set_xticks(range(0,7))
    ts.grid(linestyle=':', linewidth=0.5)
    #ts.legend(fontsize=12, loc="upper left")
    
    # # show the input surface heat fluxes for the reaction model where they are used
    # if 'BioShading_onlyC' in aqcfile:
    #     sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
    #     sx.plot(time_croco,s_flux)
    #     sx.set_xticklabels([])
    #     sx.set_ylabel('Surface heat flux (W m$^{-2}$)', fontsize=15,)
    #     sx.set_xlim(0,7)
    #     sx.set_xticks(range(0,7))
    
    plt.savefig('plot_eulerian_'+reaction+'_r'+max_photo+'_mld'+str(mld)+'_kappa'+str(kappa)+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    mlds = [20] #[10,20,50]
    kappas = [0.0001] #[0.0001,0.001,0.01,0.1]     
    reactions = ['Sverdrup_incl_K'] #['Sverdrup','Sverdrup_incl_K']
    max_photo = '10.0'
    
    for mld in mlds:
        for kappa in kappas:
            for reaction in reactions:
                do_the_plot(mld,kappa,reaction,max_photo)
    
    
