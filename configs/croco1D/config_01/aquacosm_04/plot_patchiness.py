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

def do_the_plot(mld,amplitude,mean_tau,Qswmax,React,p):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(3,1,i+1) for i in range(3)]   
    
    # get the croco input data
    crocodir='../physics/'
    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename
    time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)

    z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
    Nt_croco,Nz_croco=shape(temp_croco)
    # repeat time along the z dimension for plotting
    time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
    # repeat depth along the time dimension for plotting
    z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0)
    
    # time snapshot to plot
    t=10 #days    
            
    # plot the eulerian data
    eulfile='eulerian_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_amp'+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    time_eul,z_eul,chl_eul=get_eul_output(eulfile)
    # interpolate z_therm onto eulerian time axis (just in case different)
    z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:,0])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
    chl_eul_avg=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul)
    
    # get the aquacosm data
    aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_amp'+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    time_aqc,z_aqc,z_rank,chl_aqc,_ = get_aqc_output(aqcfile)
    
    # compute the patchiness of the aquacosms
    # (normalised standard deviation in the vicinity of each aquacosm)
    rad=2.0 # search radius (m) for computation of local stdev around each aquacosm
    stdev_chl_aqc = get_aqc_stdev(chl_aqc,z_aqc,rad)
    
    # compute the normalised realised reaction
    # get r directly from the reactions function
    r = get_aqc_reactions(time_aqc*3600*24,z_aqc,z_rank,chl_aqc,wc,React)
    
    # DO THE PLOTS        
    Nt_aqc,Nz_aqc=shape(chl_aqc)
    # repeat time along the z dimension for plotting
    time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
    
    plot_C(ax,ii+1,time_aqc,z_aqc,chl_aqc,'Aquacosms, p = '+"{0:1.0e}".format(p),t,max_chl)  
    # add the thermocline
    for n in range(0,3):
        cnt=ax[n].contour(time_croco,z_croco,temp_croco,[11],colors='w',linewidths=2.5,linestyles='dashed')
    
    # compare the chlorophyll averaged over the surface layer
    # ts=gcf().add_axes((0.4, -0.8, 0.45, 0.6))
    ts=gcf().add_axes((0., 0.55-0.55*3, 0.8, 0.5))
    ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=4, label='Eulerian')
    # ps=[1e-3,1e-7]
    ps=[1e-3,1e-4,1e-5,1e-6,1e-7,0]
    for ii,p in enumerate(ps):
        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_amp'+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
        time_aqc,z_aqc,z_rank,chl_aqc,_ = get_aqc_output(aqcfile)
        chl_aqc_avg=get_Cs_aquacosm(time_croco[:,0],z_therm_croco,time_aqc,z_aqc,chl_aqc)
        ts.plot(time_aqc,chl_aqc_avg, linewidth=2, label='Aquacosms, p = '+"{0:1.0e}".format(p))
    ts.set_ylabel('average surface Chl (mg m$^{-3}$)', fontsize=15)
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
    sx.plot(time_croco[:,0],s_flux)
    sx.set_xticklabels([])
    sx.set_ylabel('$Q_s$ (W m$^{-2}$)', fontsize=15,)
    sx.set_xlim(0,21)
    sx.set_xticks(range(0,22))
    sx.set_ylim(0,800)
    
    plt.savefig('plot_eul_aqc_2_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_amp'+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    amplitudes = [0.03] #[0, 0.01, 0.02, 0.03, 0.04]
    mlds = [10] #[10, 25]   
    mean_taus = [0] #[0, 0.05]
    Qswmaxs = [250] #[0, 250, 800]  
    
    for amplitude in amplitudes:
        for mld in mlds:
            for Qswmax in Qswmaxs:
                for mean_tau in mean_taus:
                    crocodir='../physics/'
                    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
                    crocofile=crocodir+crocofilename
                    
                    dt = 5             
                    wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
                    React = set_up_reaction(wc, dt, BioShading_onlyC,
                                        LightDecay=5.,
                                        MaxPhotoRate = 1.0, 
                                        BasalMetabolism = 0.1,
                                        Chl_C = 0.017,
                                        CrowdingMortality = 0.65,
                                        CrowdingHalfSaturation = 125,
                                        Chl_light_abs = 0.)
                    p=1e-7
                    do_the_plot(mld,amplitude,mean_tau,Qswmax,React,p)
    
    
