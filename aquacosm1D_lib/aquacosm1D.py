from pylab import *
from aquacosm1D_reactions   import set_up_reaction, NoReactions, \
    Sverdrup, Sverdrup_incl_K, SimpleBFM, BioShading, BioShading_onlyC
from aquacosm1D_diffusion   import set_up_diffusion, sort_by_depth
from aquacosm1D_transport   import set_up_transport
from aquacosm1D_watercolumn import water_column, water_column_netcdf
from aquacosm1D_utilities   import Aquacosm1D_Particles, create_particles, \
    sort_by_index, sort_by_depth, average_between, gaussian_estimate_field, \
    gaussian_average_between

seterr(divide='raise')
seterr(invalid='raise')
seterr(over='warn')
seterr(under='warn')
ion()




