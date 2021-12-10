# Simulations based on schematised diffusivity

The simulations are run at different depths, each depth representing the mixed layer depth, 
so we are isolating the response of the surface layer

`config_*/physics/make_input.py` writes netcdf files filled with schematised input values in the format needed for the aquacosm simulations

configuration descriptions
---------------------------
`config_01` - using constant diffusivities and surface radiation 

`config_02` - using diurnally varying diffusivities and surface radiation 
