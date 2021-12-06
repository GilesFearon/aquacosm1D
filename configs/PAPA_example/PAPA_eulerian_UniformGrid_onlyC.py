from pylab import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.interpolate import splev
ion()

seterr(divide='raise')
seterr(invalid='raise')
seterr(over='raise')
seterr(under='warn')


dt     = 0.75  #in seconds
Ndays  = 4*365  #length of the simulation
Nloops = int(24*60*60  *  Ndays  / dt)
Nstore = int(6*60*60 / dt) #store the particles every Nshow time steps

#-----------------------------------------------------------

Eulerian_data  = Dataset('PAPA_6h_1Y_L75_DN.nc')

#Read the diffusivity profiles
DeepestNode = 32
kappaNemoGrid = Eulerian_data.variables['difvho'][:,:DeepestNode,0,0]
zNemoGrid = Eulerian_data.variables['depthw'][:]#[:DeepestNode]
zNemoGrid_t=Eulerian_data.variables['deptht'][:]#[:DeepestNode]

#setting up the light interpolator
swr = Eulerian_data.variables['rsntds'][:,0,0]
swr *= swr>0.
times  = Eulerian_data.variables['time_counter'][:]
dtimes = (times[1]-times[0])/2
etimes = concatenate([[times[0]-dtimes],
                      times,
                      [times[-1]+dtimes, times[-1]+dtimes*2]])
eswr   = concatenate([[swr[0]], swr, [swr[-1], swr[-1]]])
swr_interpolator = interp1d(etimes, eswr, kind='linear')

#Grid for the scalars
zc = arange(0., 201., 1) - 0.5
#Grid for the diffusivity
zw = arange(0., 200., 1)

#-----------------------------------------------------------
def apply_boundary_cond(C):
    """ No Flux boundary conditions """
    C[0]  = C[1]
    C[-1] = C[-2]

#-----------------------------------------------------------
def Integrate_Chl(zc, C, Chl_C_ratio):
    """Uses trapezoids to integrate the chlorophill"""
    I = zeros_like(C)
    L = C*Chl_C_ratio
    #assume uniform L from z=0 down to the first node inside the comp.domain
    I[1] = zc[1]*L[1]
    I[2:] = 0.5*(L[2:]+L[1:-1]) * (zc[2:]-zc[1:-1])
    return cumsum(I)

def calc_rhs(kappa, swr, zw, zc, C):
        
    LightDecay = 10.
    AlphaEpsilon = 1.38e-5*(0.4/0.217)
    MaxPhotoRate = 2./(60.*60.*24.)
    BasalMetabolism = 0.16/(60.*60.*24.)
    Chl_C_ratio = 0.017
    CrowdingMortality = 1./(60.*60.*24.)
    CrowdingHalfSaturation = 12.5
    Chl_light_abs = 0.03
    
    #Diffusion terms
    D2C = ( kappa[1: ]*(C[2: ]-C[1:-1])/(zc[2: ]-zc[1:-1]) 
          - kappa[:-1]*(C[1:-1]-C[:-2])/(zc[1:-1]-zc[:-2])
          )/(zw[1:]-zw[:-1])

    #Reaction terms
    Int_cL = Chl_light_abs * Integrate_Chl(zc, C, Chl_C_ratio)
    mu = AlphaEpsilon*exp(- zc/LightDecay - Int_cL)*swr
    fE = -expm1(-Chl_C_ratio*mu/MaxPhotoRate)
    dC_gpp = fE*MaxPhotoRate*C
    dC_rsp = (BasalMetabolism*C +
              CrowdingMortality*C*C/(C+CrowdingHalfSaturation))
    
    #Adding all up
    Cdot = dC_gpp - dC_rsp
    # The end points are not used and will be fixed by b.c.
    Cdot[1:-1] += D2C
    return Cdot

    
    
#------------------------------------------------------------

class bspline(object):
    """Convenience wrapper around the splev function of scipy.

Sets up a monotone cubic B-spline approximation to the eddy
diffusivity data (z=depth, kappa=eddy diffusivity).

NOTE: this removes the first and last kappa values and replaces them
with a parabolic extrapolation of the two nearest neighboring points,
so as to enforce d(kappa)/dz = 0 at the boundaries (z[0] and z[-1]), 
which makes it easier to enforce mirror boundary conditions.

    """
    def __init__(self, z, kappa):
        kappa_min = 1.e-10 # default bdry value if extrapolation is negative
        self.degree = 3
        kappa[0] = (kappa[1]*z[2]**2 - kappa[2]*z[1]**2) / (z[2]**2 - z[1]**2)
        dz3 = z[-3] - z[-1] 
        dz2 = z[-2] - z[-1] 
        kappa[-1] = (kappa[-2]*dz3**2 - kappa[-3]*dz2**2) / (dz3**2 - dz2**2)
        if kappa[0] < kappa_min:
            kappa[0] = kappa_min
        if kappa[-1] < kappa_min:
            kappa[-1] = kappa_min            
        self.kappa = r_[kappa[4:0:-1],
                        kappa.copy(),
                        kappa[-2:-6:-1]]
        self.z = r_[2*z[0]-z[4:0:-1],
                    z.copy(),
                    2*z[-1]-z[-2:-6:-1]]
        self.t = r_[(self.z[0],)*(self.degree+1),
                    self.z[2:-2] ,
                    (self.z[-1],)*(self.degree+1)]

    def __call__(self, xx):
        return splev(xx, (self.t, self.kappa, self.degree))

    def derivative(self, xx):
        return splev(xx, (self.t, self.kappa, self.degree), der=1)

    def derivative2(self, xx):
        return splev(xx, (self.t, self.kappa, self.degree), der=2)

        
#-----------------------------------------------------------
#initial conditions
C = 0.2*(zc<100) + 1e-20

Cstore = []
LoopsPerYear = int(365*24*60*60/dt)  ### Achtung if it isn't an integer!
for loop in range(Nloops):
    #resetting time every year
    if loop % LoopsPerYear == 0:
        time = etimes[0]
        old_idx = -1
        #saving partial results
        array(Cstore).tofile("PAPA_eulerianC.onlyC.pydat")
    if loop%Nstore == 0:
        Cstore.append(C.copy())
        print("Days elapsed:", loop*dt/(24*3600.),
              'Mean C[19]:', Integrate_Chl(zc, C, 1)[19]/zc[19])
    idx = find(time <= times+dtimes)[0]
    if idx > old_idx:
        old_idx = idx
        k_interp = bspline(zNemoGrid, kappaNemoGrid[idx])
        current_kappa = k_interp(zw)
    ##########################################
    # RK2
    Cdot = calc_rhs(current_kappa,
                    swr_interpolator(time),
                    zw, zc, C)
    Caux = C + Cdot*dt/2
    apply_boundary_cond(Caux)
    Cdot = calc_rhs(current_kappa,
                    swr_interpolator(time+dt/2),
                    zw, zc, Caux)
    C += Cdot*dt
    apply_boundary_cond(C)
    ##########################################

    time += dt


Cstore = array(Cstore)
Cstore.tofile("PAPA_eulerianC.onlyC.pydat")


