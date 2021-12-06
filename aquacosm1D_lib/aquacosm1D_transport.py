from pylab import *

#------------------------------------------------------------
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
#------------------------------------------------------------

def random_walk_MilsteinStiff(z, wc, dt):
    """From 
Yin and Gan Advances in Difference Equations (2015) 2015:369
DOI 10.1186/s13662-015-0699-9
    """
    def impose_reflecting_bc(zz, wc):
        itop = find(zz < 0)
        ibtm = find(zz > wc.max_depth)
        zz[itop] = (               - zz[itop])
        zz[ibtm] = (2*wc.max_depth - zz[ibtm])

    N01   = randn(len(z))
    K     = wc.diffusivity(z)
    Kp    = wc.diffusivity_derivative(z)
    zaux  = z + N01*sqrt(2*K*dt) + 0.5*Kp*(N01*N01 + 1)*dt
    Kpaux = wc.diffusivity_derivative(zaux)
    Kppx  = wc.diffusivity_second_derivative(zaux)
    znew  = zaux + dt*(Kpaux - Kp)/(1 - Kppx)
    impose_reflecting_bc(znew, wc)

    
    return znew

#---------------------------------------------------------------------

def random_walk_Heun(z, wc, dt):
    """
    """
    def impose_reflecting_bc(zz, wc):
        itop = find(zz < 0)
        ibtm = find(zz > wc.max_depth)
        zz[itop] = (               - zz[itop])
        zz[ibtm] = (2*wc.max_depth - zz[ibtm])

    N01   = randn(len(z))
    K     = wc.diffusivity(z)
    Kp    = wc.diffusivity_derivative(z)
    b     = sqrt(2*K*dt)
    zaux  = z + Kp*dt + b*N01
    Kpaux = wc.diffusivity_derivative(zaux)
    znew  = z + 0.5*(Kp + Kpaux)*dt + b*N01 
    impose_reflecting_bc(znew, wc)
    
    return znew

#---------------------------------------------------------------------

def random_walk_Milstein(z, wc, dt):
    """
    """
    def impose_reflecting_bc(zz, wc):
        itop = find(zz < 0)
        ibtm = find(zz > wc.max_depth)
        zz[itop] = (               - zz[itop])
        zz[ibtm] = (2*wc.max_depth - zz[ibtm])

    N01 = randn(len(z))
    K  = wc.diffusivity(z)
    Kp = wc.diffusivity_derivative(z)
    znew = z + N01*sqrt(2*K*dt) + 0.5*Kp*(N01*N01 + 1)*dt
    impose_reflecting_bc(znew, wc)
    
    return znew

#---------------------------------------------------------------------

def random_walk_Visser(z, wc, dt):
    """Returns a new z as a random update of the old z after a time dt,
based on the diffusivity specified in wc. Reflecting boundary
conditions are applied.  wc must be a water_column instance.  The
algorithm follows Visser (1997) Marine Ecology Progress Series 158,
275-281.

z  : the 1D array of the particles' depths
wc : a water_column instance

Returns an updated 1D array of the particles' depths

    """
    def impose_reflecting_bc(zz, wc):
        itop = find(zz < 0)
        ibtm = find(zz > wc.max_depth)
        zz[itop] = (               - zz[itop])
        zz[ibtm] = (2*wc.max_depth - zz[ibtm])

        
    Kp = wc.diffusivity_derivative(z)
    zaux = z + dt*Kp/2.
    znew = z + dt*Kp + randn(len(zaux))*sqrt(2*wc.diffusivity(zaux)*dt)
    impose_reflecting_bc(znew, wc)
    
    return znew


#---------------------------------------------------------------------

class set_up_transport:
    """Randomly moves the particles in such a way as to mimic the eddy
diffusion profile specified in WaterColumn (a water_column instance)
"""
    def __init__(self, WaterColumn, time_step, SODE='Milstein'):
        self.wc = WaterColumn
        self.dt = time_step
        if SODE == 'Visser':
            self._sode = random_walk_Visser
        elif SODE == 'Milstein':
            self._sode = random_walk_Milstein
        elif SODE == 'Heun':
            self._sode = random_walk_Heun
        elif SODE == 'MilsteinStiff':
            self._sode = random_walk_MilsteinStiff
        #elif SODE == 'Strong1.5':
        #    self._sode = random_walk_Strong1_5
        else:
            raise ValueError("Scheme <"+SODE+"> is unknown!")
            
    def __call__(self, Particles):
        z = Particles[:,1]
        Particles[:,1] = self._sode(z, self.wc, self.dt)
        Particles.depth_sorted = False


#---------------------------------------------------------------------
