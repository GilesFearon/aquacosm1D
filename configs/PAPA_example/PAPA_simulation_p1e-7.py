import sys
sys.path.insert(0, '../aquacosm1D_lib')
from aquacosm1D import *

seed(1234567)



#------------------------------------------------------------
dt        = 5. # time step in seconds
Ndays     = 4*365 #length of the simulation
Nloops    = int(24*60*60  *  Ndays  / dt)
Nstore    = int(6*60*60 / dt) #store the particles every Nshow time steps
Npts      = 200  #number of particles
Nscalars  = 1    #number of scalars carried by each particle
Pstore    = []   #list that stores the particles every Nstore time steps

wc = water_column_netcdf(DatasetName="PAPA_6h_1Y_L75_DN.nc", max_depth=200)


Diffuse = set_up_diffusion(Npts, Nscalars,
                           radius=10.,
                           p=1.e-7, ###################
                           dt=dt,
                           wc=wc)

Transport = set_up_transport(wc, dt, SODE='Milstein')

React = set_up_reaction(wc, dt, BioShading_onlyC,
                        LightDecay=10.,
                        MaxPhotoRate = 2., 
                        BasalMetabolism = 0.16,
                        Chl_C = 0.017,
                        CrowdingMortality = 1.)

Particles = create_particles(Npts, Nscalars, wc)
zz = 100 - Particles[:,1]
Particles[:, 2] = 0.2*(zz>0) + 1e-20

starttime = wc.times[0]
wc.set_current_time(starttime)
for loop in range(Nloops):
    #-------------------------
    if loop % Nstore == 0:
        sort_by_depth(Particles)
        Pstore.append(Particles.copy())
        print('Days elapsed: ', loop*dt/(24*3600.), 'Mean C',
              mean(Particles[:,2])
        )
    Diffuse(Particles)
    React(Particles)
    Transport(Particles)

    try:
        wc.increment_current_time(dt)
    except ValueError:
        #Save partial results
        Pst = array(Pstore)
        Pst.tofile("PAPA_p{0:1.0e}.Milstein_dt5.onlyC.pydat".format(Diffuse.p))
        wc.set_current_time(starttime)

    

Pstore = array(Pstore)
Pstore.tofile("PAPA_p{0:1.0e}.Milstein_dt5.onlyC.pydat".format(Diffuse.p))


