import sys
sys.path.insert(0, '../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
ion()

def get_aqc_output(aqcfile):
    # get the aquacosm output
    data_aqc=Dataset(aqcfile)
    p=data_aqc.couping_parameter_p
    Chl_C=np.float(data_aqc.Chl_C) # ratio of chlorophyll to carbon [mgChl/mgC]
    chl_aqc=data_aqc.variables['chl'][:,:]
    z_aqc=data_aqc.variables['depth'][:,:]
    z_rank=data_aqc.variables['depth_rank'][:]
    time_aqc=data_aqc.variables['time'][:]/3600/24 # time in days
    time_aqc=time_aqc-time_aqc[0]
    chl_aqc_avg = np.mean(chl_aqc, axis=1)
    Nt_aqc,Nz_aqc=shape(chl_aqc)
    # repeat time along the z dimension for plotting
    time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
    data_aqc.close()
    return time_aqc,z_aqc,chl_aqc,chl_aqc_avg

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

def do_the_plot(mld,kappa,reaction,max_photo):
    
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(3,1,i+1) for i in range(3)]   
    
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
    
    # plot the aquacosm data
    ps=[1e-3,1e-7]
    for ii,p in enumerate(ps):
        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_'+reaction+'_r'+max_photo+'_mld'+str(mld)+"_kappa"+str(kappa)+".nc"   
        time_aqc,z_aqc,chl_aqc,chl_aqc_avg = get_aqc_output(aqcfile)
        plot_C(ax,ii+1,time_aqc,z_aqc,chl_aqc,'Aquacosms, p = '+"{0:1.0e}".format(p),t,max_chl)  
    
    # tindx_aqc = (np.abs(time_aqc[:,0] - t)).argmin()
    
    # # compare the profiles at a given time
    # bx=gcf().add_axes((0.0, -0.8, 0.3, 0.6))
    # #
    # z_eul_t = z_eul[tindx_eul,:]
    # chl_eul_t = chl_eul[tindx_eul,:]
    # bx.plot(chl_eul_t, z_eul_t, '-r', linewidth=3, label="Eulerian")
    # #
    # z_aqc_t = z_aqc[tindx_aqc,:]
    # chl_aqc_t = chl_aqc[tindx_aqc,:]
    # Particles_for_gaussian=array([z_rank,z_aqc_t,chl_aqc_t]).transpose() # reconstructing the particle format for parsing to the gaussian function
    # z_aqc_gaus, chl_aqc_gaus = gaussian_estimate_field(Particles_for_gaussian,2, mld, 1.25) # default stddev = 2.5
    # bx.plot(chl_aqc_gaus, z_aqc_gaus, '-k', linewidth=3, label="Coarse-grained")
    # imgr = bx.scatter(chl_aqc_t, z_aqc_t, c=chl_aqc_t, cmap=cm.viridis, s=15,zorder=3,label="Aquacosms")
    # imgr.set_clim(0, max_chl)
    # bx.set_ylim(51, -1)
    # bx.set_xlim(0-0.025*max_chl,max_chl+0.025*max_chl)
    # bx.set_xlabel('Chlorophyll  (mg m$^{-3}$)', fontsize=15)
    # bx.set_ylabel('Depth (m)', fontsize=15)
    # bx.text(0,-2,'Time = '+str(t)+' days',ha='left',fontsize=12)
    # bx.grid(linestyle=':', linewidth=0.5)
    # bx.legend(fontsize=12, loc="lower right") #, bbox_to_anchor=(0.1, 0.75),framealpha=0.99)
    
    # compare the chlorophyll averaged over the surface layer
    # ts=gcf().add_axes((0.4, -0.8, 0.45, 0.6))
    ts=gcf().add_axes((0., 0.55-0.55*3, 0.8, 0.5))
    ts.plot(time_eul[:,0],chl_eul_avg, 'k', linewidth=4, label='Eulerian')
    ps=[1e-3,1e-5,1e-7,1e-9]
    for ii,p in enumerate(ps):
        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_'+reaction+'_r'+max_photo+'_mld'+str(mld)+"_kappa"+str(kappa)+".nc"   
        time_aqc,z_aqc,chl_aqc,chl_aqc_avg = get_aqc_output(aqcfile)
        ts.plot(time_aqc[:,0],chl_aqc_avg, linewidth=2, label='Aquacosms, p = '+"{0:1.0e}".format(p))
    ts.set_ylabel('average surafce Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    ts.set_ylim(4, max_chl)
    ts.set_xlim(0,7)
    ts.set_xticks(range(0,7))
    ts.grid(linestyle=':', linewidth=0.5)
    ts.legend(fontsize=12, loc="upper left")
    
    # # show the input surface heat fluxes for the reaction model where they are used
    # if 'BioShading_onlyC' in aqcfile:
    #     sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
    #     sx.plot(time_croco,s_flux)
    #     sx.set_xticklabels([])
    #     sx.set_ylabel('Surface heat flux (W m$^{-2}$)', fontsize=15,)
    #     sx.set_xlim(0,7)
    #     sx.set_xticks(range(0,7))
    
    plt.savefig('plot_eul_aqc_'+reaction+'_r'+max_photo+'_mld'+str(mld)+'_kappa'+str(kappa)+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    mlds = [20,50] #[20,50]
    kappas = [0.0001] #[0.0001,0.001,0.01]     
    reactions = ['Sverdrup_incl_K'] #['Sverdrup','Sverdrup_incl_K']
    max_photo = '10.0'
    
    for mld in mlds:
        for kappa in kappas:
            for reaction in reactions:
                do_the_plot(mld,kappa,reaction,max_photo)
    
    