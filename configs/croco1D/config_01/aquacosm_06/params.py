# trying this out to simplify file naming
# (I really like this approach... could expand by including the reaction model too
# and could have a separate class which defines time step, number of aquacosms etc
# so a scalable system which would allow the testing of many different configurations
# without getting the file names too confusing while documenting what was done)
class reactions01:
    Name = 'reactions01'
    LightDecay = 5.
    AlphaEpsilon = 1.38e-5*(0.4/0.217)
    MaxPhotoRate = 1.
    BasalMetabolism = 0.16
    Chl_C = 0.017
    CrowdingMortality = 0.65
    CrowdingHalfSaturation = 125
    Chl_light_abs = 0.
    
class physics01:
    Name = 'physics01'
    amplitude = 0.03
    mean_tau = 0
    mld = 10
    Qswmax = 250
    
class physics02:
    Name = 'physics02'
    amplitude = 0.04
    mean_tau = 0
    mld = 10
    Qswmax = 250
