physics - the croco 1D simulations used to force the aquacosm runs

aquacosm_** dirs all have run_eulerian.py and run_aquacosm.py files which do what you'd expect from the name

aquacosm_01 - running for 21 days, chl ini = 1 in surface layer, tests with NoReactions, Sverdrup, Sverdrup_incl_K

aquacosm_02 - as per 01, tests with Bioshading_onlyC

aquacosm_03 - as per 02, sensitivity tests to time-step

aquacosm_04 - as per 03, CrowdingHalfSaturation = 125

aquacosm_05 - tests on multiple populations

aquacosm_06 - as per 04, setting chl ini = 0.001 in the subsurface (was 10^-20)
              used to plot summary of R depth and chl over the runs
