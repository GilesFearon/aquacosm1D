###########################################
# Imports
###########################################
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import sys
sys.path.append('/home/osboxes/2101_postdoc/aquacosm1D/croco_1D/scm_oce/')
from scm_oce import scm_oce
from scipy.io import netcdf
from scipy.interpolate import interp1d
import xarray as xr

###########################################
# Functions
###########################################
# 
########################################################
# Definition de la grille verticale du modele oceanique
# ->  z_r  : profondeur au centre des mailles
# ->  z_w  : profondeur aux interfaces
# ->  Hz   : epaisseur des mailles
########################################################
def define_oce_grid(nz,theta_s,theta_b,hc,hmax): 
    global z_r
    global z_w
    global Hz
# GGF edited function - assuming vtransform = 2 i.e. new vertical sigma coords
# for direct comparison with 2D and 3D model results. Method taken from
# matlab croco_tools function zlevs.m
#    
    z_r    = np.zeros(nz  )
    z_w    = np.zeros(nz+1)
    Hz     = np.zeros(nz)   
    z_w[0] = -hmax   
    ds  = 1./nz
    h2  = hmax+hc
    #cff = (hmax-hc)/np.sinh(theta_s)
    for k in range(nz,0,-1):
        # w grid
        sc_w     = ds * float(k-nz)
        csrf_w   = (1-np.cosh(sc_w*theta_s))/(np.cosh(theta_s)-1)
        cs_w     = (np.exp(theta_b*csrf_w)-1)/(1-np.exp(-1*theta_b))
        z0_w     = hc*sc_w + cs_w*hmax
        z_w[k  ] = z0_w*hmax/h2
        #z_w[k  ] = hc*sc_w + cff*np.sinh(theta_s*sc_w)   
        # r grid
        sc_r     = ds*(float(k-nz)-0.5)  
        csrf_r   = (1-np.cosh(sc_r*theta_s))/(np.cosh(theta_s)-1)
        cs_r     = (np.exp(theta_b*csrf_r)-1)/(1-np.exp(-1*theta_b))
        z0_r     = hc*sc_r + cs_r*hmax
        z_r[k-1] = z0_r*hmax/h2
        #z_r[k-1]   = hc*sc_r + cff*np.sinh(theta_s*sc_r)      
#
    for k in range(nz):
        Hz[k] = z_w[k+1]-z_w[k]
###########################################


##############################################################
# Definition d'un forcage par la grande echelle 
# avec un temps de rappel delta associe 
# ignore a partir du moment ou delta = 0
##############################################################
def define_largescale_oce(nz,ntra):  
    global unudge
    global vnudge    
    global tnudge     
    global delta     

    unudge = np.zeros( nz )
    vnudge = np.zeros( nz )
    tnudge = np.zeros((nz,ntra), order='F')  
    delta  = np.zeros((nz,ntra), order='F')    

##############################################################
# Definition des conditions initiales 
# uoce    : courant zonal [m/s]
# voce    : courant meridien [m/s] 
# toce[ntra=0] : temperature [Celsius]
# toce[ntra=1] : salinite [psu]
##############################################################
def define_ini_oce(nz,nt,ntra,T0,alpha,N0,mld):    
    global uoce
    global voce
    global toce

# Courant nuls    
    uoce  = np.zeros((nz,nt), order='F')
    voce  = np.zeros((nz,nt), order='F')

# Salinite constante = 35    
    toce  = np.zeros((nz,nt,ntra), order='F') 
    for t in range(nt):
        toce[:,t,1] = 35.

# Temperature  
    # Stratification ArcTH
        T1=10.
        z0=-1*mld
        deltaz0=3
        for t in range(nt):
            for k in range(nz):         
                toce[k,t,0] = 0.5*(T0-T1)*np.tanh((z_r[k]-z0)/deltaz0) + 0.5*(T0+T1)

# Tracer
    if (ntra>2):
        # ini used in Fearon et al (2020)
        # temp_mld=0.5*(T0+T1) # temperature at the defined MLD is by definition midway between T0 and T1
        # temp_target=[-999, 10, temp_mld, 999]
        # tracer_target=[1, 1, 0, 0]
        # for t in range(nt):
        #     for k in range(nz):
        #         toce[k,t,2] = interp1d(temp_target,tracer_target)(toce[k,t,0])

        # ini used for comparison with aquacosms
        for t in range(nt):
            z0=interp1d(toce[:,t,0],z_r)(11)
            for k in range(nz):         
                if z_r[k] > z0: # z is negative to > means closer to the surface than z0
                    toce[k,t,2] = 1.
                else:
                    toce[k,t,2] = 0.
        
##############################################################
# Initilisation des parametres pour la fermeture turbulente
# !!! Ne pas modifier !!! 
# nuwm   : valeur minimum pour la viscosite turbulente  
# nuws   : valeur minimum pour la diffusivite turbulente
##############################################################
def define_ini_turb(nz,ntra,ngls,nuwm,nuws,turbulence_scheme,stability_function):
    global Akt
    global Akv
    global turb
    #=====       
    tke_min = 1.e-06
    eps_min = 1.e-12
    #=====       
    Akv  = np.zeros(nz+1)  
    Akt  = np.zeros((nz+1,ntra), order='F')
    turb = np.zeros((nz+1,2,ngls), order='F')
    #=====       
    Akv [:]     = nuwm
    Akt [:,0]   = nuws
    Akt [:,1]   = nuws
    #=====       
    Akv [0]       = 0.
    Akv [nz]      = 0.
    Akt [0,0]     = 0.
    Akt [nz,0]    = 0.
    Akt [0,0]     = 0.
    Akt [nz,0]    = 0.        
    #=====   
    rp    =  0.0 ; rm    = 0.0  ; rn     =  0.0        
    if   turbulence_scheme == 1:
        rp    = -1.0 ; rm    = 0.5  ; rn     = -1.0
    elif turbulence_scheme == 2:
        rp    = 3.0  ; rm    = 1.5  ; rn     = -1.0
    elif turbulence_scheme == 3:
        rp    = 0.0  ; rm    = 1.0  ; rn     = -0.67 
    #=====              
    c1=5.;c2=0.8;c3=1.968;c4=1.136
    if   stability_function == 1:
        c1=3.6;c2=0.8;c3=1.2;c4=1.2
    elif stability_function == 2:
        c1=6.;c2=0.32;c3=0.;c4=0.
    elif stability_function == 3:
        c1=6.;c2=0.32;c3=0.;c4=0.
    elif stability_function == 4:
        c1=3.;c2=0.8;c3=2.;c4=1.118
    elif stability_function == 5:
        c1=5.;c2=0.6983;c3=1.9664;c4=1.094
    elif stability_function == 6:          
        c1=5.;c2=0.7983;c3=1.968;c4=1.136    
    #=====  
    nn  = 0.5*c1
    a1  = 0.66666666667-0.5*c2
    a2  = 1.-0.5*c3
    a3  = 1.-0.5*c4            
    cm0 =  pow( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn), 0.25 )     
    cff     = pow(cm0,3)*pow(tke_min,1.5) / eps_min                     
    gls_min = pow(cm0,rp)*pow(tke_min,rm) * pow( cff, rn )     
    #=====             
    turb[:,:,0] = tke_min
    turb[:,:,1] = gls_min    

###########################################
# Main
###########################################
if __name__ == "__main__":

    #======================================
    # Define variable lists to loop through
    #======================================
    # each combination of these variables will be modelled
    mlds = [10] #[5, 10, 15, 20, 25, 30]
    lats = [-30] #[-36, -34, -32, -30, -28, -26, -24]
    amplitudes = [0.01] #[0.01, 0.02, 0.03, 0.04]
    mean_taus = [0, 0.05] #[0, 0.05]
    T0s = [16] #[12, 14, 16, 18, 20]
    #hmaxs = [50] #[20, 50, 100, 200]
    hmax = 50
    Qswmaxs = [0, 250, 800]    
 
    for mld in mlds:
        for lat in lats:
            for amplitude in amplitudes:
                for T0 in T0s:
                    #for hmax in hmaxs:
                    for Qswmax in Qswmaxs:
                        for mean_tau in mean_taus:
                        
                            if mld > hmax:
                                continue

                            print("working on mean_tau"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_lat"+"%0.0f"%np.abs(lat)+"_T0"+str(T0)+"_Qswmax"+str(Qswmax))

                        #===============================================
                        # 1 - Parametres d entree
                        #===============================================
                            ntra         = 3       # nombre de traceurs : 2 (temperature + salinite)
                            nt           = 2       # nombre d'instants temporels stockes : 2 (n et n+1)
                            ngls         = 2       # nombre de variables pour la fermeture GLS : 2 (TKE + GLS)

                            rho0         = 1024.   # densite de l'ocean [kg/m3]
                            cp           = 3985.   # chaleur specifique [J / kg / K] 
                            #lat          = -30 # -32.291  # measurements at -32.291
                            fcor         = 4*np.pi/86400*np.sin(lat*np.pi/180.)      # parametre de Coriolis [1/s]
                            omega        = 2*np.pi/86400. # diurnal frequency for analytical forcing
                            dt           = 10. # 360.   # pas de temps du modele [s]
                            vonKar       = 0.4
                            #r_D          = 0.      # calculated below at each time-step
                            Zob          = 0.1     # Gildas recommendation
                            Cdb_min      = 0.0025  # Gildas recommendation
                            Cdb_max      = 0.02    # Gildas recommendation
                            Neu_bot      = True    # type de condition limite au fond 

                            lin_eos      = False    # type d'equation d'etatT0
                            #T0           = 18.     # temperature de surface du profil initial [C]
                            alpha        = 0.0002  # coefficient de l'équation d'etat pour la temperature [1/C]  
                            N0           = np.sqrt(0.01*alpha*9.81)  # frequence de Brunt Vaisala pour la definition de la condition initiale
                            N0=0.01
                            heatloss     = 0. #75.        #<-- refroidissement a la surface de l'ocean [W/m^2] 
                            #Qswmax       = 800. #264.91     #<-- valeur maximum du flux solaire radiatif [W/m^2]
                            nuwm         = 1.0e-4     #<-- valeur minimum pour la viscosite turbulente 
                            nuws         = 0.1e-4     #<-- valeur minimum pour la diffusivite turbulente 
                            nb_steps     = int(528. * 3600./dt) + 1       # nombre de pas de temps  192 heures = 8 jours
                            output_freq  = int( 30. *   60./dt)           # frequence de stockage de la solution (toutes les heures ici)  
                            nb_steps_6hr = int(6. * 3600./dt)         # number of timesteps in 6 hours

                        #===============================================
                        # 2 - definition de la grille verticale  [z_r,z_w,Hz]
                        #===============================================    
                            nz          = hmax     # adjust the number of levels based on the depth to avoid grid dependent results 
                            #hmax        = 60.    # profondeur totale du domaine de calcul [m] 
                            # NEW croco sigma grid implemented
                            # using identical settings to 2D and 3D models
                            hc          = 200.   # 
                            theta_s     = 7.     # 
                            theta_b     = 2.

                            define_oce_grid(nz,theta_s,theta_b,hc,hmax) #on definit la grille verticale a partir des parametres precedents 

                        #===============================================
                        # 3 - definition du forcage grande echelle et des conditions limites
                        #===============================================
                            define_largescale_oce(nz,ntra) 
                            # condition limite au fond pour traceurs    
                            dTdz_bot     = np.zeros(ntra)
                            #dTdz_bot[0]  = - N0*N0/(alpha*9.81)       #<-- derivee verticale constante = celle de la condition initiale [Celsius/m]
                            #dTdz_bot[1]  = 0.                         #<-- flux nulle pour salinite   [psu/m]

                        #===============================================
                        # 4 - boucle en temps
                        #===============================================
                            # Choix du schema de fermeture turbulente
                            # 
                            # Schema KPP       :  turbulence_scheme = 0; stability_function = 0 (parametre ignore); Ric = 0.3 (doit etre entre 0.15 et 0.45)
                            # Schema k-epsilon :  turbulence_scheme = 2; stability_function = 0; Ricstability_function = 0.3  
                            #    
                            turbulence_scheme  = 2
                            stability_function = 0
                            Ric                = 0.3

                            nout         = int ( (nb_steps-1)/output_freq  + 1 ) # nombre d'instants stockes pendant la simulation
                            temp_final   = np.zeros(nz)                     # stockage de la temperature a l'instant final
                            hbl_store    = np.zeros(nout)                   # stockage de la profondeur de couche limite
                            srflx_store  = np.zeros(nout)                   # stockage du flux solaire   
                            Akv_store    = np.zeros((nz+1,nout))            # stockage de Akv en cours de simulation
                            Akt_store    = np.zeros((nz+1,nout))            # stockage de Akv en cours de simulation
                            calendar     = np.zeros(nout)                   # stockage du temps des sorties
                            toce_store   = np.zeros((nz,nout))              # stockage de temperature en cours de simulation
                            tpas_store   = np.zeros((nz,nout))
                            uoce_store   = np.zeros((nz,nout))              # stockage de u en cours de simulation
                            voce_store   = np.zeros((nz,nout))              # stockage de v en cours de simulation
                            tau_x_store  = np.zeros(nout)
                            tau_y_store  = np.zeros(nout)               
                            tau_bot_x_store  = np.zeros(nout)
                            calendar_nb     = np.zeros(nb_steps)
                            tau_bot_y_store  = np.zeros(nb_steps) # every timestep
                            #dpdx_store   = np.zeros(nout)
                            dzetadx_simpson_store   = np.zeros(nout)
                            dzetadx_store= np.zeros(nout)

                            # Definition des conditions initiales [uoce,voce,toce]
                            define_ini_oce(nz,nt,ntra,T0,alpha,N0,mld)
                            # Initialisation des parametres du schema de fermeture [Akt,Akv,turb]
                            define_ini_turb(nz,ntra,ngls,nuwm,nuws,turbulence_scheme,stability_function)

                            # Integration temporelle 
                            rho1        = np.zeros(nz  )  # on initialise la densite a zero
                            hbls        = np.array(0.)    # on initialise la profondeur de couche limite a zero
                            ncount      = 0               # compteur pour le stockage en cours de simulation      
                            #===============================================               
                            # Boucle temporelle
                            #===============================================  
                            for kt in range(nb_steps):

                                # get time info
                                nstp = 1 +  kt    % 2     # indice temporel de l'instant n 
                                nnew = 1 + (kt+1) % 2     # indice temporel de l'instant n+1   
                                time = dt*float(kt)       # temps en secondes 
                                time_day  = time/3600/24   # time in days

                                # options for the analytical wind stress
                                #amplitude=0.02 # of rotating wind stress (N/m2)
                                anticlock=1 # 1=anticlockwise rotation, 0=clockwise rotation
                                
                                # ramp up to mean winds over a day to avoid generation of spurious oscillations 
                                if time_day < 1:
                                    mean_tau_now=time_day*mean_tau
                                else:
                                    mean_tau_now=mean_tau

                                # compute the wind stress
                                # tau_x_phase = np.pi # phase of pi means tau_x peaks at mid-day
                                tau_x_phase = 2.6619 # phase of diurnal anticlockwise rotary component of measured winds at the Nortier station over 7-14 March 2011 (see Giles /phd/data/SAWS/processed/rotary/meas_nortier_stress_rot_7day.m)
                                if anticlock==1:
                                    tau_y_phase = tau_x_phase - np.pi/2 
                                else:
                                    tau_y_phase = tau_x_phase + np.pi/2
                                tau_x_now = amplitude*np.cos(omega*time+tau_x_phase)
                                tau_y_now = amplitude*np.cos(omega*time+tau_y_phase)+mean_tau_now
                                i_tau_y_now = amplitude*np.sin(omega*time+tau_y_phase)+mean_tau_now#the sin() here effectively means we use the wind that was blowing 0.25 days ago (phase pi/2 before)    
                                                                
                                # stress input to model
                                sustr    = tau_x_now/rho0
                                svstr    = tau_y_now/rho0

                                # options for surface flux
                                #flx_flag = 1 
                                #if flx_flag==0: # zero surface fluxes
                                #    srflx    = 0.
                                #    heatflx  = 0.
                                #elif flx_flag==1: # analytical surface fluxes
                                #    srflx    = max(np.cos(omega*time - np.pi),0.)  * Qswmax / (rho0*cp)  # flux solaire (variation diurne) [C m / s]
                                #    heatflx  = srflx - heatloss / (rho0*cp)     # flux de chaleur nette  [C m / s]
                                #
                                if Qswmax == 250: # for this specific case we're assuming a constant heat flux
                                    srflx = Qswmax/(rho0*cp)
                                else:
                                    srflx    = max(np.cos(omega*time - np.pi),0.)  * Qswmax / (rho0*cp)  # flux solaire (variation diurne) [C m / s]
                                heatflx  = srflx - heatloss / (rho0*cp)     # flux de chaleur nette  [C m / s]

                                # ignoring fresh flux in these runs
                                freshflx = 0.

                                # bottom drag
                                u1= uoce[0,nstp-1]
                                v1= voce[0,nstp-1]
                                cff= np.sqrt( u1*u1 + v1*v1 )
                                # r_D = cff*pow(vonKar/((1.+Zob/Hz[0])*np.log(1.+Hz[0]/Zob) -1.),2) # from user manual
                                # croco 3D code for bottom friction:
                                Cdb = pow(vonKar/np.log(Hz[0]/Zob),2)
                                Cdb = min(Cdb_max,max(Cdb_min,Cdb))
                                r_D = cff*Cdb 
                                #r_D=0. # set friction to zero
                                tau_bot_x_now=r_D*u1
                                tau_bot_y_now=r_D*v1


                                # options for dpdx
                                dpdx_flag = 1
                                if dpdx_flag==0:
                                    dpdx_simpson=0.
                                    dpdx=dpdx_simpson
                                elif dpdx_flag==1:
                                    # note dpdx as calculated here is g*detadx
                                    dpdx_simpson=(tau_x_now+fcor/omega*i_tau_y_now)/(rho0*hmax) # Simpson (2002) condition, after Craig 1989
                                    # now add the effect of bottom friction
                                    if kt >= nb_steps_6hr:
                                        # as we don't know the solution for vb apriori we assume a diurnal frequency and use the friction 6 hours prior (pi/2 phase assuming diurnal signal)
                                        i_tau_bot_y_now=tau_bot_y_store[kt-nb_steps_6hr]*9.81 # need to mult by g here because we divided by g when we saved it
                                    else:
                                        i_tau_bot_y_now=0
                                    dpdx = dpdx_simpson-tau_bot_x_now/hmax - fcor/omega*i_tau_bot_y_now/hmax 
                                elif dpdx_flag==2:
                                    # dpdx as per dpdx_simpson above (no bottom friction effect)
                                    dpdx_simpson=(tau_x_now+fcor/omega*i_tau_y_now)/(rho0*hmax)
                                    dpdx=dpdx_simpson

                                # integration temporelle                                            
                                scm_oce.obl_stp(z_r,z_w,Hz,
                                                unudge,vnudge,tnudge,
                                                uoce,voce,toce,turb,
                                                rho0,rho1,Akv,Akt,r_D,
                                                sustr,svstr,srflx,heatflx,freshflx,
                                                dTdz_bot,delta,fcor,Ric,hbls,
                                                dt,dpdx,turbulence_scheme,stability_function,
                                                lin_eos,alpha,T0,Zob,Neu_bot,
                                                nstp,nnew,nz,ntra,nt,ngls)

                                tau_bot_y_store[kt]      = tau_bot_y_now / 9.81 # saving every timestep
                                calendar_nb[kt]          = time / (24. * 3600.)
                                # si on est a un instant de sortie
                                if  kt % output_freq == 0:
                                    calendar [ncount]        = time / (24. * 3600.)   # temps courant          
                                    hbl_store[ncount]        = -1 * hbls              # stockage de hbls
                                    srflx_store[ncount]      = heatflx * (rho0*cp)      # stockage du flux solaire  
                                    Akv_store[:,ncount]      = Akv[:]                 # stockage de Akv
                                    Akt_store[:,ncount]      = Akt[:,0]                 # stockage de Akv
                                    toce_store[:,ncount]      = toce[:,nnew-1,0]
                                    tpas_store[:,ncount]      = toce[:,nnew-1,2]
                                    uoce_store[:,ncount]      = uoce[:,nnew-1]
                                    voce_store[:,ncount]      = voce[:,nnew-1]
                                    tau_x_store[ncount]      = tau_x_now
                                    tau_y_store[ncount]      = tau_y_now
                                    tau_bot_x_store[ncount]      = tau_bot_x_now / 9.81
                                    dzetadx_simpson_store[ncount] = dpdx_simpson / 9.81
                                    dzetadx_store[ncount]       = dpdx / 9.81
                                    ncount = ncount + 1      
                            #===============================================  
                            temp_final[:] = toce[:,nnew-1,0]  # stockage du profil de temperature final    
                            del uoce
                            del voce
                            del toce
                            del turb
                            del Akv
                            del Akt            

                            # save output to .mat file

#                            scipy.io.savemat("tur"+str(turbulence_scheme)+"_mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+
#                                "_acw"+str(anticlock)+"_dzeta"+str(dpdx_flag)+"_flx"+str(flx_flag)+"_lat"+"%0.0f"%np.abs(lat)+"_T0"+str(T0)+"_hmax"+str(hmax)+".mat",
#                                    {'tau_x': tau_x_store, 'tau_y': tau_y_store, 'dzeta_dx': dzetadx_store, 'dzeta_dx_simpson': dzetadx_simpson_store,'u': uoce_store,'v': voce_store,'tpas': tpas_store,
#                                        'temp': toce_store, 'srflx': srflx_store, 'z_r': z_r, 'z_w': z_w, 'time': calendar, 'akv': Akv_store, 'akt': Akt_store})
                            
                            # write a netcdf in the same format at the NEMO netcdf file which the aquacosm model is hard coded to ingest
                            # start by reshuffling the arrays
                            calendar=calendar*24.*3600. # now in seconds
                            z_w=-np.flipud(z_w) # depth increases positively from surface downward
                            z_r=-np.flipud(z_r) # depth increases positively from surface downward
                            Akt_store=np.flipud(Akt_store)
                            toce_store=np.flipud(toce_store)
                            tpas_store=np.flipud(tpas_store)
                            uoce_store=np.flipud(uoce_store)
                            voce_store=np.flipud(voce_store)
                            Akt_store=np.expand_dims(np.expand_dims(Akt_store.transpose(),axis=-1),axis=-1)
                            tpas_store=np.expand_dims(np.expand_dims(tpas_store.transpose(),axis=-1),axis=-1)
                            toce_store=np.expand_dims(np.expand_dims(toce_store.transpose(),axis=-1),axis=-1)
                            uoce_store=np.expand_dims(np.expand_dims(uoce_store.transpose(),axis=-1),axis=-1)
                            voce_store=np.expand_dims(np.expand_dims(voce_store.transpose(),axis=-1),axis=-1)
                            srflx_store=np.expand_dims(np.expand_dims(srflx_store,axis=-1),axis=-1)
                            tau_x_store=np.expand_dims(np.expand_dims(tau_x_store,axis=-1),axis=-1)
                            tau_y_store=np.expand_dims(np.expand_dims(tau_y_store,axis=-1),axis=-1)
                            dzetadx_store=np.expand_dims(np.expand_dims(dzetadx_store,axis=-1),axis=-1)
                            ds = xr.Dataset(data_vars=dict(difvho=(["time_counter","depthw","y","x"],Akt_store),
                                                                   temp=(["time_counter","deptht","y","x"],toce_store),
                                                                   tpas=(["time_counter","deptht","y","x"],tpas_store),
                                                                   u=(["time_counter","deptht","y","x"],uoce_store),
                                                                   v=(["time_counter","deptht","y","x"],voce_store),
                                                                   tau_x=(["time_counter","y","x"],tau_x_store),
                                                                   tau_y=(["time_counter","y","x"],tau_y_store),
                                                                   dzetadx=(["time_counter","y","x"],dzetadx_store),
                                                                   rsntds=(["time_counter","y","x"],srflx_store),
                                                                   ),
                                                                   coords=dict(depthw=(["depthw"],z_w),deptht=(["deptht"],z_r),time_counter=(["time_counter"],calendar),x=(["x"],[0]),y=(["y"],[0])))
                            ds.time_counter.attrs['units'] = 'seconds since 0000-1-1 00:00:00' # dummy time as an analytical experiment
                            # add global attributes telling us what the input parameters were
                            ds.attrs["tau_ac"] = str(amplitude)+' N m^-2'
                            ds.attrs["tau_mean"] = str(mean_tau)+' N m^-2'
                            ds.attrs["mld"] = str(mld)+' m'
                            ds.attrs["lat"] = str(lat)+' deg'
                            ds.attrs["ini_surface_temp"] = str(T0)+' deg C'
                            ds.attrs["depth"] = str(hmax)+' m'
                            # write the output
                            fname_out="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat"+"%0.0f"%np.abs(lat)+"_T0"+str(T0)+"_hmax"+str(hmax)+".nc"
                            ds.to_netcdf(fname_out)

    print("Done!")


