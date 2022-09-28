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
    if n==2:
        # ax[n].set_xlabel('Time (days)', fontsize=15,)
        cbarax = gcf().add_axes((0.82, -0.5, 0.015, 1.5)) # [left, bottom, width, height] 
        cbar = colorbar(img, cbarax)
        cbar.set_label("Chlorophyll  (mg m$^{-3}$)", fontsize=15)
    # else:
    #     ax[n].set_xticklabels([])
    
    ax[n].set_xticklabels([])
    
    ax[n].text(0.2, 45, label,
               fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    # ax[n].plot([t, t], [0, 50], '--k',
    #                     linewidth=1)

def do_the_plot(mld,kappa,React):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(3,1,i+1) for i in range(3)]   
    
    # get the croco input data
    physicsdir='../physics/'
    physicsfilename="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
    physicsfile=physicsdir+physicsfilename
    time_physics, z, zw, kappa_physics, s_flux = get_physics_input(physicsfile)
    
    # time snapshot to plot
    t=10 #days    
            
    # plot the eulerian data
    eulfile='eulerian_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mld'+str(mld)+'_kappa'+str(kappa)+'_dt10.0.nc'
    time_eul,z_eul,chl_eul,chl_eul_avg=get_eul_output(eulfile)
    # repeat time along the z dimension for plotting
    Nt_eul,Nz_eul=shape(chl_eul)
    time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # repeat depth along the time dimension for plotting
    z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    tindx_eul = (np.abs(time_eul[:,0] - t)).argmin()
    #
    # max values to plot
    max_chl = 1.5
    plot_C(ax,0,time_eul,z_eul,chl_eul,'Eulerian',t,max_chl)
    
    # plot the aquacosm data
    if React.dt == 5.:
        ps=[1e-3,1e-7]
    else:
        # ps=[2e-3,2e-7]
        ps=[2e-4,2e-8]
    for ii,p in enumerate(ps): 
        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mld'+str(mld)+'_kappa'+str(kappa)+'_dt'+str(React.dt)+'.nc'
        time_aqc,z_aqc,z_rank,chl_aqc,_ = get_aqc_output(aqcfile)
        Nt_aqc,Nz_aqc=shape(chl_aqc)
        # repeat time along the z dimension for plotting
        time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
        plot_C(ax,ii+1,time_aqc,z_aqc,chl_aqc,'Aquacosms, p = '+"{0:1.0e}".format(p),t,max_chl)  
    
    # compare the chlorophyll averaged over the surface layer
    # ts=gcf().add_axes((0.4, -0.8, 0.45, 0.6))
    ts=gcf().add_axes((0., 0.55-0.55*3, 0.8, 0.5))
    ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=4, label='Eulerian')
    if React.dt == 5.:
        ps=[1e-3,1e-7]
    else:
        # ps=[2e-3,2e-7]
        ps=[2e-4,2e-8]
    for ii,p in enumerate(ps):
        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mld'+str(mld)+'_kappa'+str(kappa)+'_dt'+str(React.dt)+'.nc'
        time_aqc,z_aqc,z_rank,chl_aqc,chl_aqc_avg = get_aqc_output(aqcfile)
        ts.plot(time_aqc,chl_aqc_avg, linewidth=2, label='Aquacosms, p = '+"{0:1.0e}".format(p))
    ts.set_ylabel('average Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    # max values to plot
    # base = 2 #nearest multiple of 'base'
    # max_chl = base * round(np.max(chl_eul_avg)/base)
    # max_chl = 12#15#35#500
    ts.set_ylim(0,max_chl)
    ts.set_xlim(0,21)
    ts.set_xticks(range(0,22))
    ts.grid(linestyle=':', linewidth=0.5)
    ts.legend(fontsize=12, loc="upper left")
    
    # # show the input surface radiation
    sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
    sx.plot(time_physics,s_flux)
    sx.set_xticklabels([])
    sx.set_ylabel('$Q_s$ (W m$^{-2}$)', fontsize=15,)
    sx.set_xlim(0,21)
    sx.set_xticks(range(0,22))
    sx.set_ylim(0,800)
    
    plt.savefig('plot_eul_aqc_2_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mld'+str(mld)+'_kappa'+str(kappa)+'_dt'+str(React.dt)+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":

    mlds = [20]  
    kappas = [0.0001] #[0.0001,0.001,0.01]  
    
    for kappa in kappas:
        for mld in mlds:
            crocodir='../physics/'
            crocofilename="mld"+str(mld)+"_kappa"+str(kappa)+".nc"
            crocofile=crocodir+crocofilename
            
            dt = 5.        
            wc = water_column_netcdf(DatasetName=crocofile, max_depth=mld)
            React = set_up_reaction(wc, dt, BioShading_onlyC,
                                LightDecay=5.,
                                MaxPhotoRate = 1., 
                                BasalMetabolism = 0.1,
                                Chl_C = 0.017,
                                CrowdingMortality = 0.25,
                                Chl_light_abs = 0.)
            
            do_the_plot(mld,kappa,React)

    
