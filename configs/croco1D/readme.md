# Simulations based on 1D CROCO model

`scm_oce` contains the Fortran source code for the 1D croco model

`docs` has some description of the model (see also Fearon et al., 2020)

We can call this code in python by first running f2py to create a python library from the Fortran code
You need to first create and activate a python 3.4 environment- needed for f2py to work.
See https://conda.io/docs/user-guide/tasks/manage-environments.html
I've called mine py34, so 'conda activate py34' to activate the environment
Then 'cd scm_oce' and compile the f2py library- do this by 'make clean' then 'make'
This creates a python library from the fortran 1D source code

`config_*/physics/run_croco.py` runs the 1D code (do this in your py34 environment) and writes netcdf files in the format needed for the aquacosm simulations

configuration descriptions
---------------------------
`config_01` - analytically configured model of Fearon et al. (2020) i.e. landsea breeze forcing at 30degS, but initialised with a passive tracer with a value of 1 in surface and zero in the subsurface, and run duration extended to 21 days

`config_02` - (placerholder for other croco1D setups) 
