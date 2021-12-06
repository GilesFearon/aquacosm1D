import sys
sys.path.insert(0, '../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
ion()


DatasetName = 'PAPA_6h_1Y_L75_DN.nc'
Eul  = Dataset(DatasetName)

#--------------------------------------------------------------------
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
#--------------------------------------------------------------------

Ntimes = 4*365*4
Npts = 200
Nfields = 3
SavesPerDay = 4


zl  = 1
Crb = 2
LoC = 0.017

days = arange(Ntimes, dtype='float')/SavesPerDay
ddays = datetime(2011, 1, 1, 0, 0) + days*timedelta(days=1)

wc = water_column_netcdf(DatasetName="PAPA_6h_1Y_L75_DN.nc", max_depth=200)
start_time = wc.get_current_time()
#-----------------------------------------------------------
zcU = arange(0., 201., 1) - 0.5
MyeulUnifC = reshape(fromfile('PAPA_eulerianC.onlyC.pydat'),
                     (Ntimes, 201))

#---

# ParticleFiles = [
#     "PAPA_p1e-02.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-03.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-04.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-05.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-06.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-07.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-08.Milstein_dt5.onlyC.pydat",
#     "PAPA_p1e-09.Milstein_dt5.onlyC.pydat",
#     "PAPA_p0e+00.Milstein_dt5.onlyC.pydat",
#     ]

ParticleFiles = [
    "PAPA_p1e-07.Milstein_dt5.onlyC.pydat",
    ]

# Labels = ["p=1e-2", "p=1.e-3", "p=1.e-4", "p=1.e-5", "p=1.e-6",
#           "p=1.e-7", "p=1.e-8", "p=1.e-9", "p=0"]

Labels = ["p=1.e-7"]

Pstores = []
for Fname in ParticleFiles:
    Pstores.append(reshape(fromfile(Fname), (Ntimes, Npts, Nfields)))


#-------------------------------------------------
P0    = Pstores[-1]
P1em7 = Pstores[-4]
P1em3 = Pstores[ 1]
YearOffset = 4*365*3
Nframes = YearOffset+array([440, 608, 800, 1032])
Pcontent = {
    "Eulerian": MyeulUnifC,
    "p=1e-3": P1em3,
    "p=1e-7": P1em7,
    "p=0": P0,
}
figure(figsize=(9,10))
ax = [subplot(4,1,i+1) for i in range(4)]

for n, label in enumerate(Pcontent):
    ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
    times = []
    depths = []
    conc = []
    if label == "Eulerian":
        for i, Nodes in enumerate(Pcontent[label]):
            times.append(ddays[i]+ timedelta(days=1)*ones((len(zcU),)) )
            depths.append(zcU)
            conc.append(Nodes*LoC)
    else:
        for i, PP in enumerate(Pcontent[label]):
            times.append(ddays[i]+ timedelta(days=1)*ones((Npts,)) )
            depths.append(PP[:,1])
            conc.append(PP[:,Crb]*LoC)
    img = ax[n].scatter(times, depths, c=log10(conc),
                        cmap=cm.viridis, linewidth=0, s=40)
    img.set_clim(-5, 1)
    if label == "Eulerian":
        mld = Eul.variables["mldr10_1"][:,0,0]
        ax[n].plot(ddays[-len(mld):], mld, '-w', linewidth=0.5)
    ax[n].set_ylim(200, 0)
    ax[n].set_xlim(datetime(2014,1,1,0,0),
                   datetime(2015,1,1,0,0))
    ax[n].text(datetime(2014,1,4,0,0), 185, label,
               fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    ax[n].set_xticks([datetime(2014, 1,1,0,0),
                      datetime(2014, 2,1,0,0),
                      datetime(2014, 3,1,0,0),
                      datetime(2014, 4,1,0,0),
                      datetime(2014, 5,1,0,0),
                      datetime(2014, 6,1,0,0),
                      datetime(2014, 7,1,0,0),
                      datetime(2014, 8,1,0,0),
                      datetime(2014, 9,1,0,0),
                      datetime(2014,10,1,0,0),
                      datetime(2014,11,1,0,0),
                      datetime(2014,12,1,0,0),
                      datetime(2015, 1,1,0,0),
    ])
    if n<3:
        ax[n].set_xticklabels([])
    else:
        ax[n].set_xticklabels(["             Jan",
                               "             Feb", "             Mar",
                               "             Apr", "             May",
                               "             Jun", "             Jul",
                               "             Aug", "             Sep",
                               "             Oct", "             Nov",
                               "             Dec", 
                               ""])
gcf().text(0.01, 0.76, "Depth (m)", fontsize=15, rotation=90)
ax[-1].set_xlabel("Time", fontsize=15)
cbarax = gcf().add_axes((0.89, 0.44, 0.02, 0.55))
cbar = colorbar(img, cbarax)
cbar.set_label("Chlorophyll  log$_{10}$(mg/m$^3$)", fontsize=15)

#--
transparentax=gcf().add_axes((0.08, 0.44, 0.8, 0.55))
transparentax.patch.set_alpha(0)
for key in transparentax.spines:
    transparentax.spines[key].set_visible(False)
transparentax.set_xticks([])
transparentax.set_yticks([])
for Nfr in Nframes:
    transparentax.plot([ddays[Nfr], ddays[Nfr]], [0, 1], '--k',
                       linewidth=1)
transparentax.set_ylim(0, 1)
transparentax.set_xlim(datetime(2014,1,1,0,0),
                       datetime(2015,1,1,0,0))
#--

bx = [[gcf().add_axes((0.08+j*0.23, 0.265-i*0.1, 0.18, 0.1))
       for j in range(4)]
      for i in range(3)]
for i in range(3):
    for j in range(4):
        if i<2:
            bx[i][j].set_xticklabels([])
        #if j>0:
        #    bx[i][j].set_yticklabels([])
gcf().text(0.01, 0.255, "Depth (m)", fontsize=15, rotation=90)
gcf().text(0.4, 0.01, "Chlorophyll (mg/m$^3$)", fontsize=15)

ylims = [130, 40, 25, 49]
xlims = [0.07, 2, 7, 0.5]
for j, Nfr in enumerate(Nframes):
    for n, label in enumerate(Pcontent):
        if label=="Eulerian":
            continue
        i = n-1
        z = Pcontent[label][Nfr][:,1]
        P = Pcontent[label][Nfr][:,Crb]*LoC
        zsm, csm = gaussian_estimate_field(Pcontent[label][Nfr],
                                           Crb, 200.0, 2.5)
        mld = Eul.variables["mldr10_1"][Nfr-YearOffset,0,0]
        wc.set_current_time(start_time+(Nfr-YearOffset)*6*3600)
        kp = wc.diffusivity_derivative(wc.z)
        zmix = wc.z[find(kp==amin(kp))[0]]
        bx[i][j].plot([-1,10], [zmix, zmix], '--k', zorder=-1)
        bx[i][j].plot([-1,10], [mld,   mld], '-', color='gray', zorder=-1)
        bx[i][j].plot(Pcontent["Eulerian"][Nfr]*LoC, zcU,
                      '-r', linewidth=3, label="Eulerian")
        bx[i][j].plot(csm*LoC, zsm, '-k', linewidth=3,
                      label="Coarse-grained")
        imgr = bx[i][j].scatter(P, z, c=P, cmap=cm.viridis, s=15,
                                #zorder=3,
                                label="Aquacosms")
        bx[i][j].set_ylim(ylims[j], -1)
        bx[i][j].set_xlim(-0.05*xlims[j], xlims[j])
        imgr.set_clim(0, xlims[j]*0.8)
        bx[i][j].text(xlims[j]*0.4, ylims[j]*0.97, label, fontsize=9)
bx[0][3].legend(fontsize=7, loc="center",
                bbox_to_anchor=(0.1, 0.75),
                framealpha=0.99)
#--
Nframes2 = YearOffset+array([165, 590, 1015, 1442])
tr2ax=gcf().add_axes((0.08, 0.365, 0.8, 0.075))
tr2ax.patch.set_alpha(0)
for key in tr2ax.spines:
    tr2ax.spines[key].set_visible(False)
tr2ax.set_xticks([])
tr2ax.set_yticks([])
for Nfr, Nfr2 in zip(Nframes, Nframes2):
    tr2ax.plot([ddays[Nfr], ddays[Nfr2]], [0, 1], '--k',
                       linewidth=1)
tr2ax.set_ylim(1, 0)
tr2ax.set_xlim(datetime(2014,1,1,0,0),
                       datetime(2015,1,1,0,0))
