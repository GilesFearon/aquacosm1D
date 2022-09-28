import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
ion()
    
def do_the_anim(mld,amplitude,mean_tau,Qswmax,React,p):
    
    fname='aquacosm_p'+"{0:1.0e}".format(p)+'_'+type(React.current_model).__name__+'_'+'mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)
    aqcfile=fname+'.nc'
    print('\n working on ' + aqcfile +'\n')
        
    # get the croco input data
    crocodir='../physics/'
    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename
    time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
    #
    z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco,11)
    Nt_croco,Nz_croco=shape(temp_croco)
    # repeat time along the z dimension for plotting
    time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
    # repeat depth along the time dimension for plotting
    z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0) 
            
    # get the eulerian data
    eulfile='eulerian_'+type(React.current_model).__name__+'_'+'mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    
    time_eul,z_eul,chl_eul=get_eul_output(eulfile)
    # interpolate z_therm onto eulerian time axis (just in case different)
    z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:,0])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
    chl_eul_avg=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul)
    #
    # repeat time along the z dimension for plotting
    Nt_eul,Nz_eul=shape(chl_eul)
    time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # repeat depth along the time dimension for plotting
    z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    
    # get the aquacosm data
    # aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_amp'+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    time_aqc,z_aqc,z_rank,chl_aqc,_ = get_aqc_output(aqcfile)
    Nt_aqc,Nz_aqc=shape(chl_aqc)
    # get average chl in surface layer for aquacosm run
    chl_aqc_avg=get_Cs_aquacosm(time_croco[:,0],z_therm_croco,time_aqc,z_aqc,chl_aqc)
    # repeat time along the z dimension for plotting
    time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)

    # get the diagnostic data 
    diagfile=fname+'_diags.nc'
    _,_,_,_,r,_,_,_,_=get_aqc_diags(diagfile) # only need normalised growth rate r
    
    # do the plot
    fig=figure(figsize=(10,12))
    #
    # max values to plot
    max_chl = 12
    #
    # plot the aquacosm data
    ax1 = subplot(1,1,1) 
    ax1.set_position(  (0.08, 0.75, 0.8, 0.23))
    #
    # plot_C(ax1,time_aqc,z_aqc,chl_aqc,max_chl)  
    img = ax1.scatter(time_aqc, z_aqc, c=chl_aqc,
                        cmap=cm.viridis, linewidth=0, s=40)
    img.set_clim(0, max_chl)
    ax1.set_ylim(50, 0)
    ax1.set_xlim(0,21)
    ax1.set_xticks(range(0,22))
    ax1.set_ylabel('Depth (m)', fontsize=15,)
    ax1.set_xticklabels([])   
    # add the colorar
    cbarax = gcf().add_axes((0.9, 0.75, 0.015, 0.23)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    cbar.set_label("Chlorophyll  (mg m$^{-3}$)", fontsize=15)
    # add the thermocline
    cnt=ax1.contour(time_croco,z_croco,temp_croco,[11],colors='w',linewidths=2.5,linestyles='dashed')
    
    # time-series plot
    ts=gcf().add_axes((0.08, 0.49, 0.8, 0.23))
    ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=2, label='Eulerian')
    ts.plot(time_aqc[:,0],chl_aqc_avg, linewidth=2, label='Aquacosms, p = '+"{0:1.0e}".format(p))
    ts.set_ylabel('average surface Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    ts.set_ylim(0,max_chl)
    ts.set_xlim(0,21)
    ts.set_xticks(range(0,22))
    ts.grid(linestyle=':', linewidth=0.5)
    ts.legend(fontsize=12, loc="upper left")
    
    # time snapshot to plot
    t=0 #days  
    tindx_aqc = (np.abs(time_aqc[:,0] - t)).argmin()
    tindx_eul = (np.abs(time_eul[:,0] - t)).argmin()
    
    # compare the profiles at a given time
    bx=gcf().add_axes((0.2, 0.05, 0.55, 0.37))
    #
    z_eul_t = z_eul[tindx_eul,:]
    chl_eul_t = chl_eul[tindx_eul,:]
    eul_line, = bx.plot(chl_eul_t, z_eul_t, '-r', linewidth=3, label="Eulerian")
    #
    z_aqc_t = z_aqc[tindx_aqc,:]
    chl_aqc_t = chl_aqc[tindx_aqc,:]
    Particles_for_gaussian=array([z_rank,z_aqc_t,chl_aqc_t]).transpose() # reconstructing the particle format for parsing to the gaussian function
    z_aqc_gaus, chl_aqc_gaus = gaussian_estimate_field(Particles_for_gaussian,2, 50, 1.25) # default stddev = 2.5
    aqc_coarse, = bx.plot(chl_aqc_gaus, z_aqc_gaus, '-k', linewidth=3, label="Coarse-grained")
    aqc_scat = bx.scatter(chl_aqc_t, z_aqc_t, c=chl_aqc_t, cmap=cm.viridis, s=15,zorder=3,label="Aquacosms")
    
    aqc_scat.set_clim(0, max_chl)
    bx.set_ylim(51, -1)
    max_chl_xlim = 30
    bx.set_xlim(0-0.025*max_chl_xlim,max_chl_xlim+0.025*max_chl_xlim)
    bx.set_xlabel('Chlorophyll  (mg m$^{-3}$)', fontsize=15)
    bx.set_ylabel('Depth (m)', fontsize=15)
    bx_title=bx.text(max_chl_xlim/2,-2,'Time = '+str(t)+' days',ha='center',fontsize=12)
    bx.grid(linestyle=':', linewidth=0.5)
    bx.legend(fontsize=12, loc="lower right") #, bbox_to_anchor=(0.1, 0.75),framealpha=0.99)
    
    ts_timestamp,=ts.plot([t,t],[-999,999],color='orange',linestyle=':', linewidth=2)
    ax1_timestamp,=ax1.plot([t,t],[-999,999],color='orange',linestyle=':',  linewidth=2)
    
    def animate(i):
        tindx_aqc = (np.abs(time_aqc[:,0] - i)).argmin()
        tindx_eul = (np.abs(time_eul[:,0] - i)).argmin()
        z_eul_t = z_eul[tindx_eul,:]
        chl_eul_t = chl_eul[tindx_eul,:]
        z_aqc_t = z_aqc[tindx_aqc,:]
        chl_aqc_t = chl_aqc[tindx_aqc,:]
        Particles_for_gaussian=array([z_rank,z_aqc_t,chl_aqc_t]).transpose() # reconstructing the particle format for parsing to the gaussian function
        z_aqc_gaus, chl_aqc_gaus = gaussian_estimate_field(Particles_for_gaussian,2, 50, 1.25) # default stddev = 2.5
        
        eul_line.set_xdata(chl_eul_t)
        
        aqc_coarse.set_ydata(z_aqc_gaus)
        aqc_coarse.set_xdata(chl_aqc_gaus)
        
        aqc_scat.set_offsets(np.c_[chl_aqc_t, z_aqc_t])
        aqc_scat.set_array(chl_aqc_t)
        
        ts_timestamp.set_xdata([i,i])
        ax1_timestamp.set_xdata([i,i])
        
        # bx_title.set_text('Time = '+str(i)+' days')
        bx_title.set_text('Time = {0:.2f} days'.format(i))
    
    anim = FuncAnimation(fig, animate, frames=np.arange(0,21,6/24))
    
    # really struggled to install imagemagick on the WSL Ubunutu system I'm working on so using PillowWriter instead
    # seems the only down side is that it doesn't automatically set the extent of the saved image like "bbox_inches = 'tight'" when saving a jpeg
    # so you have to play around with setting the size of the figure upfront and placing the axes within the defined figure extents
    my_writer=PillowWriter(fps=10) 
    anim.save(fname+'.gif', writer=my_writer)
    # anim.save(fname+'.gif', writer='imagemagick')
    
if __name__ == "__main__":
    
    amplitudes = [0.03] #[0, 0.01, 0.02, 0.03, 0.04]
    mlds = [10] #[10, 25]
    mean_taus = [0] #[0, 0.05]
    Qswmaxs = [250] #[0, 250, 800]
    ps=[1e-7]
    for amplitude in amplitudes:
        for mld in mlds:
            for Qswmax in Qswmaxs:
                for mean_tau in mean_taus:
                    crocodir='../physics/'
                    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
                    crocofile=crocodir+crocofilename
                    
                    dt = 5             
                    wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
                    React = set_up_reaction(wc, dt, Sverdrup_incl_K, 
                                            LightDecay = 5.,
                                            BasePhotoRate = 1.,
                                            RespirationRate = 0.1,
                                            CarryingCapacity = 20)
                    
                    for p in ps:
                        do_the_anim(mld,amplitude,mean_tau,Qswmax,React,p)
    
    
