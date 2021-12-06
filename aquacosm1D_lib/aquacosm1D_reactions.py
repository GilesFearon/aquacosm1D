from pylab import *
from aquacosm1D_utilities import sort_by_depth

"""This module contains the reaction terms for the 'aquacosm1D'
program. The usage follows this pattern:
1) initialize the reaction terms as in:
React = set_up_reaction(<model name>, <put parameters here>)
2) within the main loop pass the instance 'React' to the numerical 
   ODE integrator (right now the integrator is RK2).

The models implemented so far are:

- NoReactions : scalars carried by the particles are left unchanged 
                  upon calling the ODE integrator.

- Sverdrup : Sverdrup's classic 1953 model.

- SimpleBFM : A very simple version of the Biogeochemical Flux Model tracking
                only organic carbon and chlorophyll.

- BioShading : Same as SimpleBFM but accounts for the light absorption of the 
               plankton present in the water column.

"""


#--------------------------------------------------------------

def RK2(Particles, Reaction, dt, wc):
    """Integrates the reaction terms with a second-order Runge-Kutta scheme.
    'Reaction' should be an instance of the reaction classes.
    'Particles' must have as many scalar fields as 'Reaction' returns.
    'dt' is the time step.
    'wc' is an instance of 'water_column'.

    This modifies the scalars carried by 'Particles' in-place.
    """
    aux  = Particles.copy()
    time = wc.get_current_time()
    aux[:,2:]       += Reaction(Particles, wc, time)*dt/2.
    Particles[:,2:] += Reaction(aux, wc, time+dt/2.)*dt

#--------------------------------------------------------------

class set_up_reaction:
    """Initializes the reaction terms.

It must be called as:

React = set_up_reaction(<water_column instance>,
                        <time step>,
                        <name of reaction model>, <parameters of that model>)

Then the reaction rates may be obtained by calling the 'React' instance as:

React.current_model(Particles, React.wc)

where: 

'Particles' is the array of particles' depths and scalar concentrations
'React.wc' is the 'water_column' instance stored in React at initialization: 
           note that the time at which anything time-dependent (e.g. surface 
           radiation) is evaluated is React.wc.get_current_time()

"""
    def __init__(self, WaterColumn, dt, model, **prms):
        self.current_model = model(**prms)
        self.wc = WaterColumn
        self.dt = dt
        try:
            self.__dict__.update(self.current_model.__dict__)
        except:
            raise ValueError('Model <'+repr(model)+'> not implemented yet!')

    def __call__(self, Particles):
        """Advances by one time step the concentration fields in
'Particles'. Right now it uses a Runge-Kutta II order scheme.

        """
        RK2(Particles, self.current_model, self.dt, self.wc)

    
#------------------------------------------------------------------------

class NoReactions:
    def __init__(self, **prms):
        pass

    def __call__(self, Particles, wc, time):
        RRates = 0*Particles[:,2:]
        return RRates


#-------------------------------------------------------------------------

class Sverdrup:
    """
Sverdrups' 1953 model. See also Kida & Ito, JGR Oceans 122, 9160 (2017).
There is only one scalar P (Chlorophyll/Phytoplankton concentration).
The ODEs are:

dP/dt = (BasePhotoRate*exp(-z/LightDecay) - RespirationRate) * P

where z is the depth.

Parameters:
  LightDecay    [m]      : light absorption e-folding scale.
  BasePhotoRate [1/day]  : maximum photosynthetic rate.
  RespirationRate [1/day]: name says it all. 
"""
    def __init__(self, LightDecay = 5.,
                       BasePhotoRate = 1.,
                       RespirationRate = 0.1):
        self.LightDecay      = LightDecay
        self.BasePhotoRate   = BasePhotoRate/(60.*60.*24.)   #[1/seconds]
        self.RespirationRate = RespirationRate/(60.*60.*24.) #[1/seconds]

    def __call__(self, Particles, wc, time):
        """Uses Particles[:,2] as the chlorophyll field and Particles[:,1] as
        depth.
        """
        mu = self.BasePhotoRate*exp(-Particles[:,1]/self.LightDecay)
        RRates = 0*Particles[:,2:]
        RRates[:,0] = (mu - self.RespirationRate)*Particles[:,2]
        
        # GGF debugging - to be deleted
        check = (mu - self.RespirationRate)
        
        return RRates

#-------------------------------------------------------------------------

class Sverdrup_incl_K:
    """
As per class Sverdrup, but including a carrying capacity to limit growth.
The ODEs are:

dP/dt = (BasePhotoRate*exp(-z/LightDecay) - RespirationRate) * P * (1 - P/CarryingCapacity)

where z is the depth.

Parameters:
  LightDecay    [m]      : light absorption e-folding scale.
  BasePhotoRate [1/day]  : maximum photosynthetic rate.
  RespirationRate [1/day]: name says it all. 
  CarryingCapacity [concentration] : Maximum Chlorophyll/Phytoplankton concentration
"""
    def __init__(self, LightDecay = 5.,
                       BasePhotoRate = 1.,
                       RespirationRate = 0.1,
                       CarryingCapacity = 20):
        self.LightDecay       = LightDecay
        self.BasePhotoRate    = BasePhotoRate/(60.*60.*24.)   #[1/seconds]
        self.RespirationRate  = RespirationRate/(60.*60.*24.) #[1/seconds]
        self.CarryingCapacity = CarryingCapacity

    def __call__(self, Particles, wc, time):
        """Uses Particles[:,2] as the chlorophyll field and Particles[:,1] as
        depth.
        """
        mu = self.BasePhotoRate*exp(-Particles[:,1]/self.LightDecay)
        RRates = 0*Particles[:,2:]
        RRates[:,0] = (mu - self.RespirationRate)*Particles[:,2]*(1 - Particles[:,2]/self.CarryingCapacity)
        return RRates



#-------------------------------------------------------------------------

class SimpleBFM:
    """Marcello's Biogeochemical Flux Model with just organic carbon (C)
and chlorophyll (L). 

Parameters:
  LightDecay    [m]        : light absorption e-folding scale.
  AlphaEpsilon [mgC/(mgChl*s)] : alpha^0_chl * epsilon_par (see notes).
  MaxPhotoRate [1/day]     : maximum photosynthetic rate.
  BasalMetabolism [1/day]  : same as respiration rate in Sverdrup's model.
  Max_Chl_C [mgChl/mgC]    : Maximum chl:C ratio.
  CrowdingMortality [1/day]: Mortality rate in the crowding term.
  CrowdingHalfSaturation [mgC]: name says it all.

    """
    def __init__(self, LightDecay = 10.,
                       AlphaEpsilon = 1.38e-5*(0.4/0.217),
                       MaxPhotoRate = 2.,
                       BasalMetabolism = 0.16,
                       Max_Chl_C = 0.025,
                       CrowdingMortality = 1.,
                       CrowdingHalfSaturation = 12.5,
    ):
        self.LightDecay      = LightDecay #[m]
        self.AlphaEpsilon    = AlphaEpsilon #[mgC/(mgChl*s)]
        self.MaxPhotoRate    = MaxPhotoRate/(60.*60.*24.)   #[1/s]
        self.BasalMetabolism = BasalMetabolism/(60.*60.*24.) #[1/s]
        self.Max_Chl_C       = Max_Chl_C #[mgChl/mgC]
        self.CrowdingMortality = CrowdingMortality/(60.*60.*24.) #[1/s]
        self.CrowdingHalfSaturation = CrowdingHalfSaturation #[mgC]

    def __call__(self, Particles, wc, time):
        """Uses Particles[:,2] as carbon (C), Particles[:,3] as chlorophyll
        (L), and Particles[:,1] as depth (z).

        """
        z = Particles[:,1]
        C = Particles[:,2]
        L = Particles[:,3]
        RRates = 0*Particles[:,2:]
        Q = wc.surface_swr(time)
        mu = self.AlphaEpsilon*exp(-z/self.LightDecay)*Q
        fE = -expm1(-(mu*L)/(self.MaxPhotoRate*C))
        dC_gpp = fE*self.MaxPhotoRate*C
        dC_rsp = (self.BasalMetabolism*C +
                  self.CrowdingMortality*C*C/(C+self.CrowdingHalfSaturation))
        rhochl = ((self.Max_Chl_C*dC_gpp) / (mu*L)  if Q>0.
                  else dC_gpp*0.)
        RRates[:,0] = dC_gpp - dC_rsp
        RRates[:,1] = rhochl*dC_gpp - dC_rsp*L/C
        
        return RRates


#-------------------------------------------------------------------------

class BioShading:
    """Modified version of 'SimpleBFM' to include changes of transparency
in the water column produced by the plankton itself. There are only
two fields: organic carbon (C) and chlorophyll (L).  

Parameters:
  LightDecay    [m]        : light absorption e-folding scale.
  AlphaEpsilon [mgC/(mgChl*s)]: alpha^0_chl * epsilon_par (see notes).
  MaxPhotoRate [1/day]     : maximum photosynthetic rate.
  BasalMetabolism [1/day]  : same as respiration rate in Sverdrup's model.
  Max_Chl_C [mgChl/mgC]    : Maximum chl:C ratio.
  CrowdingMortality [1/day]: Mortality rate in the crowding term.
  CrowdingHalfSaturation [mgC]: name says it all.
  Chl_light_abs [m^2/mgChl]: Light absorption coefficient of chlorophill.

    """
    def __init__(self, LightDecay = 10.,
                       AlphaEpsilon = 1.38e-5*(0.4/0.217),
                       MaxPhotoRate = 2.,
                       BasalMetabolism = 0.16,
                       Max_Chl_C = 0.025,
                       CrowdingMortality = 1.,
                       CrowdingHalfSaturation = 12.5,
                       Chl_light_abs = 0.03,
    ):
        self.LightDecay             = LightDecay #[m]
        self.AlphaEpsilon           = AlphaEpsilon #[mgC/(mgChl*s)]
        self.MaxPhotoRate           = MaxPhotoRate/(60.*60.*24.)   #[1/s]
        self.BasalMetabolism        = BasalMetabolism/(60.*60.*24.) #[1/s]
        self.Max_Chl_C              = Max_Chl_C #[mgChl/mgC]
        self.CrowdingMortality      = CrowdingMortality/(60.*60.*24.) #[1/s]
        self.CrowdingHalfSaturation = CrowdingHalfSaturation #[mgC]
        self.Chl_light_abs          = Chl_light_abs #[m^2/mgChl]


    def __integrate_chl__(self, Particles):
        """Integrates vertically the chlorophyll field. Returns the values of
the integral at the position of the particles. Uses the trapezoid method. Particles are assumed to be sorted by height."""
        z = Particles[:,1]
        L = Particles[:,3]
        I = zeros_like(L)
        I[0]  = L[0]*z[0] #assume uniform L from z=0 down to the first particle
        I[1:] = 0.5*(L[1:]+L[:-1]) * (z[1:]-z[:-1])
        return cumsum(I)
        
    def __call__(self, Particles, wc, time):
        """Uses Particles[:,2] as carbon (C), Particles[:,3] as chlorophyll
        (L), and Particles[:,1] as depth (z).

        """
        if not Particles.depth_sorted:
            sort_by_depth(Particles)
        z = Particles[:,1]
        C = Particles[:,2]
        L = Particles[:,3]
        RRates = 0*Particles[:,2:]
        Q = wc.surface_swr(time)
        Int_cL = self.Chl_light_abs * self.__integrate_chl__(Particles)
        mu = self.AlphaEpsilon*exp(-z/self.LightDecay - Int_cL)*Q
        fE = -expm1(-(mu*L)/(self.MaxPhotoRate*C))
        dC_gpp = fE*self.MaxPhotoRate*C
        dC_rsp = (self.BasalMetabolism*C +
                  self.CrowdingMortality*C*C/(C+self.CrowdingHalfSaturation))
        rhochl = ((self.Max_Chl_C*dC_gpp) / (mu*L)  if Q>0
                  else dC_gpp*0.)
        RRates[:,0] = dC_gpp - dC_rsp
        RRates[:,1] = rhochl*dC_gpp - dC_rsp*L/C
        
        return RRates


#-------------------------------------------------------------------------

class BioShading_onlyC:
    """Modified version of 'BioShading' which uses a constant L/C ratio.
There is only one field: organic carbon (C). Chlorophyll (L) may be
later inferred by multiplying C by the specified L/C ratio.

Parameters:
  LightDecay    [m]        : light absorption e-folding scale.
  AlphaEpsilon [mgC/(mgChl*s)]: alpha^0_chl * epsilon_par (see notes).
  MaxPhotoRate [1/day]     : maximum photosynthetic rate.
  BasalMetabolism [1/day]  : same as respiration rate in Sverdrup's model.
  Chl_C [mgChl/mgC]        : Fixed Chl:C ratio.
  CrowdingMortality [1/day]: Mortality rate in the crowding term.
  CrowdingHalfSaturation [mgC]: name says it all.
  Chl_light_abs [m^2/mgChl]: Light absorption coefficient of chlorophill.

    """
    def __init__(self, LightDecay = 10.,
                       AlphaEpsilon = 1.38e-5*(0.4/0.217),
                       MaxPhotoRate = 2.,
                       BasalMetabolism = 0.16,
                       Chl_C = 0.017,
                       CrowdingMortality = 1.,
                       CrowdingHalfSaturation = 12.5,
                       Chl_light_abs = 0.03,
    ):
        self.LightDecay             = LightDecay #[m]
        self.AlphaEpsilon           = AlphaEpsilon #[mgC/(mgChl*s)]
        self.MaxPhotoRate           = MaxPhotoRate/(60.*60.*24.)   #[1/s]
        self.BasalMetabolism        = BasalMetabolism/(60.*60.*24.) #[1/s]
        self.Chl_C                  = Chl_C #[mgChl/mgC]
        self.CrowdingMortality      = CrowdingMortality/(60.*60.*24.) #[1/s]
        self.CrowdingHalfSaturation = CrowdingHalfSaturation #[mgC]
        self.Chl_light_abs          = Chl_light_abs #[m^2/mgChl]


    def __integrate_chl__(self, Particles):
        """Integrates vertically the chlorophyll field, computed as the carbon
field times the fixed L/C ratio. Returns the values of the integral at
the position of the particles. Uses the trapezoid method. Particles
are assumed to be sorted by height.

        """
        z = Particles[:,1]
        L = Particles[:,2]*self.Chl_C
        I = zeros_like(L)
        I[0]  = L[0]*z[0] #assume uniform L from z=0 down to the first particle
        I[1:] = 0.5*(L[1:]+L[:-1]) * (z[1:]-z[:-1])
        return cumsum(I)
        
    def __call__(self, Particles, wc, time):
        """Uses Particles[:,2] as carbon (C) and Particles[:,1] as depth (z).

        """
        if not Particles.depth_sorted:
            sort_by_depth(Particles)
        z = Particles[:,1]
        C = Particles[:,2]
        RRates = 0*Particles[:,2:]
        Q = wc.surface_swr(time)
        Int_cL = self.Chl_light_abs * self.__integrate_chl__(Particles)
        mu = self.AlphaEpsilon*exp(-z/self.LightDecay - Int_cL)*Q
        fE = -expm1(-self.Chl_C*mu/self.MaxPhotoRate)
        dC_gpp = fE*self.MaxPhotoRate*C
        dC_rsp = (self.BasalMetabolism*C +
                  self.CrowdingMortality*C*C/(C+self.CrowdingHalfSaturation))
        RRates[:,0] = dC_gpp - dC_rsp
        
        return RRates
