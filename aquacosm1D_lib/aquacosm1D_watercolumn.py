from pylab import *
from netCDF4 import Dataset
from scipy.interpolate import PchipInterpolator #monotonic cubic interpolator
from scipy.interpolate import interp1d
from scipy.interpolate import splev


#------------------------------------------------------------
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
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

        
    
#------------------------------------------------------------
class water_column:
    """Instances of water_column contain eddy diffusivity and other
physical data of relevance for the water column (e.g. light at the surface).
They also keep track of time.

kappa     : (m^2/s) either a positive number or a function of depth and time.
max_depth : (m) positive number; min depth is assumed to be 0.
short_wave_radiation : (w/m^2) non-negative number or a function of time.
start_time: (s) beginning of simulation.
last_time : (s) if specified, an error will be raised if current time>last_time

"""
    def __init__(self, kappa, max_depth,
                 short_wave_radiation=0,
                 start_time=0,
                 last_time=None):
        self.set_current_time(float(start_time))
        self.start_time = float(start_time)
        if last_time is None:
            self.last_time = None
        else:
            self.last_time = float(last_time)
        self._init_depth(max_depth)
        self._init_kappa(kappa)
        self._init_swr(short_wave_radiation)
    
    def _init_depth(self, max_depth):
        if isinstance(max_depth, (int, long, float)):
            self.max_depth = float(max_depth)
        else:
            raise ValueError("'max_depth' must be a number")

    def _init_kappa(self, kappa):
        if (isinstance(kappa, (int, long, float))
            and float(kappa) > 0): 
            self._kappa_value = float(kappa)
            self.diffusivity = self._constant_diffusivity
            self.diffusivity_derivative = self._constant_diffusivity_derivative
            self.diffusivity_second_derivative = \
                self._constant_diffusivity_second_derivative
        elif callable(kappa):
            try:
                float(kappa(0, self._ctime))
            except:
                raise ValueError("'kappa' must be a function of two numbers returning a single number")
            self.diffusivity = lambda z: kappa(z, self._ctime)
            self.diffusivity_derivative = self._first_finite_differences_4th
            self.diffusivity_second_derivative = \
                self._second_finite_differences_4th
            self._diff_dz = self.max_depth*1.e-8
        else:
            raise ValueError("'kappa' must be either a number larger than zero, or a function")
        
    def _init_swr(self, short_wave_radiation):
        if (isinstance(short_wave_radiation, (int, long, float))
            and float(short_wave_radiation) >= 0):
            self._swr_value = float(short_wave_radiation)
            self.surface_swr = self._constant_swr
        elif callable(short_wave_radiation):
            try:
                float(short_wave_radiation(self._ctime))
            except:
                raise ValueError("'short_wave_radiation' must be a function of a single number returning a single number")
            self.surface_swr = short_wave_radiation
        else:
            raise ValueError("'short_wave_radiation' must be either a non-negative number or a function")
            
    def _first_finite_differences_4th(self, z):
        dd = self.diffusivity
        dz = self._diff_dz
        return (dd(z-2*dz) - 8*dd(z-dz) + 8*dd(z+dz) - dd(z+2*dz))/(12*dz)
    
    def _second_finite_differences_4th(self, z):
        dd = self.diffusivity
        dz = self._diff_dz
        return (-(dd(z-2*dz)+dd(z+2*dz)) + 16*(dd(z-dz)+dd(z+dz)) - 30*dd(z))/(12*dz*dz)
    
    def _constant_swr(self, time):
        return self._swr_value

    def _constant_diffusivity(self, z):
        return self._kappa_value * ones_like(z)
    
    def _constant_diffusivity_derivative(self, z):
        return zeros_like(z)

    def _constant_diffusivity_second_derivative(self, z):
        return zeros_like(z)

    def set_current_time(self, t):
        """Sets the current time (self._ctime) to t. The preferred unit of
measurement is seconds. Consider using self.increment_current_time()
to increment the current time by one timestep.

        """
        self._ctime = t
        
    def get_current_time(self):
        """Returns the current time"""
        return self._ctime
        
    def increment_current_time(self, dt):
        """Increments the current time from self._ctime to self.ctime+dt.
Raises an error if the new time is larger than self.last_time.

        """
        if (self.last_time is None) or (self._ctime + dt <= self.last_time):
            self.set_current_time(self._ctime + dt)
        else:
            raise ValueError ("Attempting to increment time beyond self.last_time = ", self.last_time)

            


#------------------------------------------------------------
class water_column_netcdf(water_column):
    """Instances of water_column contain eddy diffusivity and other
physical data of relevance for the water column (e.g. light at the surface).
They also keep track of time.

self.z :         z-coordinates of the Eulerian grid where the diffusivity 
                 is available.

self.times :     1D array specifying the times (in sec. elapsed since 
                 1900-01-01) at which diffusivities and other data information 
                 change available.

self.max_depth : maximum depth of the Lagrangian water
                 column. Particles won't go any deeper (reflecting
                 boundary condition) than this depth.

self.data :      contains all the data of the netcdf file specified in 
                 'DatasetName' as read by netCDF4.

self.kappa :     vertical eddy diffusivity in m^2/s, there's one profile for
                 each time entry, and one level in the profile for each z entry.

self._Kintrp :   a set of monotonic cubic interpolators in space for the data 
                 in self.kappa, one interpolator for each distinct time. 
                 Evaluate this by calling self.diffusivity(z) for any z between
                 zero and self.max_depth

self._Kderiv :   same as self._Kintrp, but for the vertical derivatives. 
                 Evaluate this by calling self.diffusivity_derivative(z) for 
                 any z between zero and self.max_depth

self._Kderiv2 :  same as self._Kintrp, but for the second vertical derivatives. 
                 Evaluate this by calling 
                 self.diffusivity_second_derivative(z) for any z between zero 
                 and self.max_depth

    """
    def __init__(self,
                 DatasetName,
                 max_depth):
        self.data = Dataset(DatasetName)

        # z is depth in m, positive down
        self.z = self.data.variables['depthw'][:]
        self.max_depth = max_depth
        if max_depth is not None:
            imax = len(find(self.z<=max_depth))
            if imax == len(self.z):
                self.z[-1] = self.max_depth
            elif ((self.max_depth - self.z[imax-1]) <
                  (self.z[imax] - self.max_depth)):
                self.z = self.z[:imax]
                self.z[-1] = self.max_depth
            else:
                self.z = self.z[:imax+1]
                self.z[-1] = self.max_depth
                
        # Times in 'time_counter' are centered with respect to
        # their time intervals, instead values in self.times mark the 
        # beginning of the time interval
        self.times  = self.data.variables['time_counter'][:]
        halfperiod  = (self.times[1]-self.times[0])/2.
        self.times -= halfperiod
        self.set_current_time(self.times[0])

        # diffusivities
        self.kappa    = self.data.variables['difvho'][:,:imax,0,0]
        self._Kintrp  = [bspline(self.z, k) for k in self.kappa] # here k is kappa through depth at each timestep, the loop creates a bspline object at each timestep
        self._Kderiv  = [k.derivative  for k in self._Kintrp]
        self._Kderiv2 = [k.derivative2 for k in self._Kintrp]

        # surface incoming shortwave radiation
        self.swr = self.data.variables['rsntds'][:,0,0]
        # The following insures that the (incoming) surface radiation
        # is always non-negative. Some datasets report tiny negative
        # numbers at night.
        self.swr *= self.swr>0.
        # Extending radiation data 1/2 time interval at the beginning and end
        etimes = concatenate([[self.times[0]],
                              self.data.variables['time_counter'][:],
                              [self.times[-1]+2*halfperiod]])
        eswr = concatenate([[self.swr[0]], self.swr, [self.swr[-1]]])
        self._swr_intrp = interp1d(etimes, eswr, kind='linear')
        self.last_time = etimes[-1]

    def set_current_time(self, t):
        """Sets the current time (._ctime) to t, in seconds elapsed since
        1900-01-01 00:00:00. Useful values range between self.times[0]
        and self.times[-1]. Sets the current discrete time (._itime)
        as the index of the latest time in the past for which there is
        diffusivity data available. Diffusivities and their
        derivatives are computed at this time.
        There is no bounds check! Consider using self.increment_current_time()

        """
        self._ctime = t
        self._itime = find(self.times<=t)[-1]
            
    def increment_current_time(self, dt):
        """Increments the current time from self._ctime to self.ctime+dt. 
Raises an error if the new time is closer than  dt from self.last_time."""
        if self._ctime + dt <= self.last_time - dt:
            self.set_current_time(self._ctime + dt)
        else:
            raise ValueError ("Attempting to increment time beyond self.last_time - dt")
            
    def diffusivity(self, z):
        """Returns the eddy diffusivity at heigth z in m^2/s. 

        """
        return self._Kintrp[self._itime](z)

    def diffusivity_derivative(self, z):
        """Returns the vertical gradient of eddy diffusivity at height z in
        m/s. 

        """
        return self._Kderiv[self._itime](z)

    def diffusivity_second_derivative(self, z):
        """Returns the second derivative in the vertical direction of eddy
        diffusivity at height z in 1/s. 

        """
        return self._Kderiv2[self._itime](z)

    def surface_swr(self, time=None):
        """Returns the surface net shortwave radiation in W/m^2 interpolated
        at time 'time', or at the current time if 'time==None'.

        """
        return (self._swr_intrp(self._ctime)
                if time is None
                else self._swr_intrp(time))

    

    
