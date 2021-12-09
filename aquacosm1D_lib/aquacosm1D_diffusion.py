from pylab import *
import ctypes
from aquacosm1D_utilities import sort_by_depth
import os

#D  = ctypes.cdll.LoadLibrary("./aqc1D_diffusion.so")
#D  = ctypes.cdll.LoadLibrary("/home/franz/Desktop/Water_column/aquacosm1D_8/aqc1D_diffusion.so")
#D  = ctypes.cdll.LoadLibrary("/home/osboxes/2101_postdoc/aquacosm1D/aquacosm1D_lib/aqc1D_diffusion.so")
dir_lib = os.path.dirname(os.path.abspath(__file__))
D  = ctypes.cdll.LoadLibrary(dir_lib+"/aqc1D_diffusion.so")

#------------------------------------------------------------

class set_up_diffusion:
    """Diffuses the scalars carried by the particles using 'method 2' of
Paparella-Popolizio J. Comp. Phys., 2018.

The strength of interaction between a particle pair is zero if the
distance of the pair is larger than 'radius', otherwise it is given by
a normalized Gaussian function multiplied by p. The standard deviation
of the Gaussian is proportional to the smallest eddy diffusivity
evaluated at the particles' position. The eddy diffusivity is evaluated
through the water_column object wc.

Given a time step 'dt', the corresponding nominal diffusivity of the
gaussian is stored in 'self.NominalDiffusivity'. An estimate of the
effective diffusivity among the particles is stored in
'self.EffectiveDiffusivity'.

    """
    def __init__(self, Npts, Nscalars, p, radius, dt, wc):
        self.Npts     = Npts
        self.Nscalars = Nscalars        
        self.p = p
        self.radius = radius
        self.dt = dt
        self.wc = wc
        #This is just an auxiliary array for diffusion calculations
        self.SumExchangedMasses = empty((Npts, Nscalars), dtype='float')
    
    def __call__(self, P):
        """This is equation (17) in Paparella-Popolizio J. Comp. Phys., 2018.
Changes in-place the particle list P: it is sorted by depth
(particle depth is assumed to be the second column in P), and the
scalars are diffused. The same diffusion is used for all scalars.

        """
        if not P.depth_sorted:
            sort_by_depth(P)
        EddyDiff = self.wc.diffusivity(P[:,1])
        D.diffuse_varD(P.ctypes.data_as(ctypes.c_void_p),
                self.SumExchangedMasses.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(self.Npts), ctypes.c_int(self.Nscalars),
                ctypes.c_double(self.p),
                ctypes.c_double(self.radius),
                ctypes.c_double(self.dt),
                EddyDiff.ctypes.data_as(ctypes.c_void_p),
        )
        
