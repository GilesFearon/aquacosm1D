# Simulations based on 1D CROCO model

`scm_oce` contains the Fortran source code for the 1D croco model
(also on Gitlab: git@gitlab.uct.ac.za:croco-uct/croco-uct.git)

We can call this code in python by first running f2py to create a python library from the Fortran code
You need to first create and activate a python 3.4 environment- needed for f2py to work.
See https://conda.io/docs/user-guide/tasks/manage-environments.html
I've called mine py34, so 'conda activate py34' to activate the environment
Then 'cd scm_oce' and compile the f2py library- do this by 'make clean' then 'make'
This creates a python library from the fortran 1D source code
Then navigate to your configuration directory (still in your py34 environment)

each configuration directory `config_*` contains a `physics` dir (where the croco 1D output is generated) and an `aquacosm` dir (where the aquacosm model is run using the croco physics as input)

configuration descriptions
---------------------------
`config_01` - analytically configured model of Fearon et al. (2020) i.e. landsea breeze forcing at 30degS, but initialised with a passive tracer with a value of 1 in surface and zero in the subsurface, and run duration extended to 21 days

`config_02` - (placerholder for other croco1D setups) 
