These case files are from the 2021 cavity flame holder paper:

Sitaraman, H., Yellapantula, S., de Frahan, M. T. H., Perry, B., 
Rood, J., Grout, R., & Day, M. (2021). Adaptive mesh based combustion simulations 
of direct fuel injection effects in a supersonic cavity flame-holder. 
Combustion and Flame, 232, 111531.

The runs from the paper are quite expensive. The baseline grid 
is 512 x 128 x 32 and simulations are performed with 3 levels (max_level=2).

When running on coarser grids, there may be numerical stability issues
depending on how the velocity field is initialized. Three options for velocity
initialization are available; setting `prob.init_type = 1` or `2` may help
simulations stably get through the initial transient. These
options apply reduced velocities in the cavity region; in the default (0)
there is a uniform initial x velocity everywhere.

First a non-reacting simulation is performed and then reacting cases for 
C1 and C2 are restarted from the non-reacting checkpoint file with 
inject_fuel=1 and appropriate injection location set by centx. The 
fuel autoignites because of the high stagnation temperature near the injection
location.

Further information is also available in the youtube video from 
APS gallery of fluid motion:
https://www.youtube.com/watch?v=QIaLfaethIw
