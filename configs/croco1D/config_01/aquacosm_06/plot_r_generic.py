import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
import xarray as xr
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import params 

ion()

def do_the_plot(wc,react_params):
    
    React = set_up_reaction(wc, dt, BioShading_onlyC,
                LightDecay=react_params.LightDecay,
                AlphaEpsilon=react_params.AlphaEpsilon,
                MaxPhotoRate = react_params.MaxPhotoRate, 
                BasalMetabolism = react_params.BasalMetabolism,
                Chl_C = react_params.Chl_C,
                CrowdingMortality = react_params.CrowdingMortality,
                CrowdingHalfSaturation = react_params.CrowdingHalfSaturation,
                Chl_light_abs = react_params.Chl_light_abs)
        
    # set up aquacosm array in correct format for calling reactions library
    depth=int(wc.max_depth)
    Npts=depth+1 # 1m increments
    Nscalars=1
    Particles = create_particles(Npts, Nscalars, React.wc)
    Particles      = Aquacosm1D_Particles(
        zeros((Npts, Nscalars+2), dtype='float64')
        )
    Particles[:,1] = np.linspace(0, depth,num=Npts)
    sort_by_depth(Particles)
    Particles[:,0] = arange(Npts)
    
    # define Chlorophyll values to compute reaction rates
    # we define Chl bins based on a log scale, then create uniformly distributed
    # values with each Chl bin, and then compute the median value of each bin. 
    # This is to recreate the method used for processing the aquacosm model output
    # where the data are also binned and the median is computed (see plot_aqc_r_summary.py)
    # Chls = np.array([0.001,0.01,0.1,1,10,100])
    Chl_bin_centres = np.array([0.001,0.01,0.1,1,10,100])
    Chl_bin_edges = np.append(Chl_bin_centres/3,Chl_bin_centres[-1]*3) # roughly in the middle of the ticks when plotted on log scale
    Nchl = len(Chl_bin_centres)
    
    R_summary = np.zeros((Npts,Nchl))
    for chl_indx in range(Nchl):
        # create 100 Chl values spanning this bin
        Chl_bin_values = np.linspace(Chl_bin_edges[chl_indx], Chl_bin_edges[chl_indx+1],num=100)
        # compute r for each Chl value over all depths 
        R_bin = np.zeros((Npts,len(Chl_bin_values)))
        for ii,Chl in enumerate(Chl_bin_values):
            # update particles array with Chl
            Chl_array=np.zeros((Npts,))+Chl
            C=Chl_array/React.Chl_C 
            Particles[:,2]=C
            RRates=React.current_model(Particles, React.wc, 0)
            R_bin[:,ii]=RRates[:,0]/Particles[:,2]
        # get the median over Chl
        R_summary[:,chl_indx] = np.median(R_bin,axis=1) 
    
    R_summary = R_summary*3600*24 # now days^-1
    
    z_bin_centres=Particles[:,1]
    z_bin_edges = np.append(z_bin_centres-0.5,z_bin_centres[-1]+0.5)
        
    # do_the_plot
        
    figure(figsize=(5,5))
    ax = subplot(1,1,1) 
    
    ax.set_position(  (0., 0., 1., 0.5))
    img = ax.pcolormesh(Chl_bin_edges, z_bin_edges, R_summary, cmap=cm.RdBu_r)
    img.set_clim(-1, 1)
    ax.set_ylim(20.5, -0.5)
    ax.set_yticks(np.linspace(0,20,num=11))

    # ax.set_xlim(0.003,300)
    ax.set_xscale('log')
    #Let's get rid of the scientific notation on the x-axis
    fmt = mpl.ticker.StrMethodFormatter("{x:g}")
    ax.xaxis.set_major_formatter(fmt)
    ax.minorticks_off()
    
    for C in Chl_bin_edges:
        ax.plot([C,C],[z_bin_edges[0],z_bin_edges[-1]],color='k',linewidth=1)
    
    ax.set_ylabel('Depth (m)', fontsize=12,)
    ax.set_xlabel('Chl mg m$^{-3}$)', fontsize=12,)
    
    cbarax = gcf().add_axes((1.02, 0., 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    cbar.set_label('R (days$^{-1}$)', fontsize=12)
    
    fname_out='plot_r_generic_'+react_params.Name+'_flx'+str(React.wc.surface_swr(0))+'.jpg'
    plt.savefig(fname_out,dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    # dummy values needed here to create a water column object 
    # only shortwave radiation is used in the calculation
    dt = 5             
    kappa = 0.001
    max_depth=20
    wc = water_column(kappa, max_depth, short_wave_radiation=250)
    
    react_params=params.reactions01 # select which reactions params to plot
    
    do_the_plot(wc,react_params)
    
