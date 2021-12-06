from pylab import *

"""This module contains the eulerian diffusion code, 
useful for comparing against the aquacosm output

"""

#--------------------------------------------------------------

def RK2(kappa,zw,zc,C,dt):
    """Integrates the tracer terms with a second-order Runge-Kutta scheme.
    """

    Cdot = eulerian_diffusion(kappa,zw,zc,C)
    Caux = C + Cdot*dt/2
    apply_boundary_cond(Caux)
    
    Cdot = eulerian_diffusion(kappa,zw,zc,Caux)
    C += Cdot*dt
    apply_boundary_cond(C)
    
    return C
    
def eulerian_diffusion(kappa,zw,zc,C):
    
    D2C = ( kappa[1: ]*(C[2: ]-C[1:-1])/(zc[2: ]-zc[1:-1]) 
          - kappa[:-1]*(C[1:-1]-C[:-2])/(zc[1:-1]-zc[:-2])
          )/(zw[1:]-zw[:-1])
    
    Cdot = 0*C
    # The end points are not used and will be fixed by b.c.
    Cdot[1:-1] += D2C
    
    return Cdot

def apply_boundary_cond(C):
    """ No Flux boundary conditions """
    C[0]  = C[1]
    C[-1] = C[-2]
    
class set_up_eulerian_diffusion:
    """Diffuses the tracer(s) on the eulerian grid.

"""
    def __init__(self, WaterColumn, dt, Nscalars):
        self.wc = WaterColumn
        self.dt = dt
        self.Nscalars = Nscalars
        zw = self.wc.z
        self.zw = zw
        dz = zw[1:]-zw[0:-1] # thickness of each layer
        zc = zw[0:-1]+0.5*dz # zc (tracer depths) defined midway between zw depths
        zc = np.concatenate(([zc[0]-dz[0]],zc,[zc[-1]+dz[-1]])) # add end points - the end points are not used and will be fixed by b.c.
        self.zc = zc

    def __call__(self, Tracers):
        """Advances by one time step the concentration fields in
'Tracers' using a Runge-Kutta II order scheme.

        """
        
        for t in range(self.Nscalars):

            Tracers[:,t+2] = RK2(self.wc.kappa[self.wc._itime,:],self.zw,self.zc,Tracers[:,t+2],self.dt)
            
            
