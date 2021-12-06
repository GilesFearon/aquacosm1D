from pylab import *



class Aquacosm1D_Particles(ndarray):
    """A subclassed ndarray with attributes that allow to track whether or not the particles are sorted.
The particles array has two indexes: rows represent particles, columns are defined as follows:

P[:,0]  - the is the particle's ID: an integer (represented as float) uniquely 
          identifying each particle.
P[:,1]  - the particle's depth.
P[:,2:] - the scalars carried by the particle (if any).

Attributes:

.depth_sorted - if 'True' the particles are sorted by depth (shallowest first)
.pid_sorted   - if 'True' the particles are sorted by PID (smallest first)

See:  https://docs.scipy.org/doc/numpy-1.15.0/user/basics.subclassing.html
"""

    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.depth_sorted = False
        obj.pid_sorted   = False
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.depth_sorted = getattr(obj, 'depth_sorted', None)
        self.pid_sorted   = getattr(obj, 'pid_sorted', None)

        
#-----------------------------------------------------------------------

def create_particles(Npts, Nscalars, wc):
    """Convenience function to create a particle array respectuful of the
aquacosm1D conventions.

Npts     : number of particles
Nscalars : number of scalars carried by each particle.
wc       : a water column instance.

Returns an (Npts, Nscalars+2) array. The first column is the particle
ID (should never be changed). The second column is the particle's
depth. Particles are sorted by their depth. All other columns are
initialized to zero.

    """
    Particles      = Aquacosm1D_Particles(
        zeros((Npts, Nscalars+2), dtype='float64')
        )
    Particles[:,1] = rand(Npts)*wc.max_depth
    sort_by_depth(Particles)
    Particles[:,0] = arange(Npts)
    Particles.depth_sorted = True
    Particles.pid_sorted = True
    return Particles

#------------------------------------------------------------

def average_between(depth1, depth2, Pstore, fieldindex):
    if depth1 > depth2:
        depth1, depth2 = depth2, depth1
    out = zeros(len(Pstore), dtype='float64')
    for n, PP in enumerate(Pstore):
        idx = find((PP[:,1]>=depth1) *
                   (PP[:,1]<=depth2))
        if len(idx)>0:
            out[n] = mean(PP[idx,fieldindex])
        else: #assign the value of the closest particle
            d1 = fabs(PP[:,1]-depth1)
            d2 = fabs(PP[:,1]-depth2)
            if amin(d1) < amin(d2):
                out[n] = PP[argmin(d1),fieldindex]
            else:
                out[n] = PP[argmin(d2),fieldindex]           
    return out
        

#------------------------------------------------------------

def gaussian_estimate_field(Particles, fieldindex, max_depth, stdev):
    Nlevels = len(Particles)*4
    dzlvl = max_depth/Nlevels
    zlevels = linspace(dzlvl/2, max_depth-dzlvl/2, Nlevels)
    estimate = zeros_like(zlevels)
    normaliz = zeros_like(zlevels)
    gdenom = 2*stdev**2
    for P in Particles:
        ee = exp(- (zlevels-P[1])**2 / gdenom)
        estimate += ee*P[fieldindex]
        normaliz += ee
    return zlevels, estimate/normaliz

#------------------------------------------------------------

def gaussian_average_between(depth1, depth2, Pstore, fieldindex, wc,
                             stdev=None):
    """Returns an estimate of the average value of the scalar indexed by
'fieldindex' between depth1 and depth2.

    """
    if depth1 > depth2:
        depth1, depth2 = depth2, depth1
    if stdev is None:
        stdev = wc.max_depth/len(Pstore[0])
    out = zeros(len(Pstore), dtype='float64')
    for n, PP in enumerate(Pstore):
        zz, qq = gaussian_estimate_field(PP, fieldindex, wc.max_depth, stdev)
        idx = find((zz>=depth1) * (zz<=depth2))
        if idx is []:
            d = fabs(zz-(depth1+depth2)/2)
            idx = [find(d==amin(d))[0]]
        out[n] = mean(qq[idx])
    return out
        

#------------------------------------------------------------
def is_sorted(x):
    """Returns true if the 1D numerical array x is sorted in ascending order"""
    min = amin(x[1:]-x[:-1])
    return min < 0


def sort_by_index(particles):
    """Sort in increasing PID (particle id) the list of particles, while
maintaining each particle's scalars. It assumes that depth is the
value in the second column.

    """
    particles[:] = particles[particles[:,0].argsort(), :]
    particles.pid_sorted = True
    if is_sorted(particles[:,1]):
        particles.depth_sorted = True
    else:
        particles.depth_sorted = False
        
    
    
def sort_by_depth(particles):
    """Sort by increasing depth the list of particles in place. Each
particle retains its PID and all the scalars that it is carying. 
Assumes that depth is the value in the second column.

    """
    particles[:] = particles[particles[:,1].argsort(), :]
    particles.depth_sorted = True
    if is_sorted(particles[:,0]):
        particles.pid_sorted = True
    else:
        particles.pid_sorted = False
        
#------------------------------------------------------------
