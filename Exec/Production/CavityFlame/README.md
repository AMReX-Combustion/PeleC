These case files are from the 2021 cavity flame holder paper:

Sitaraman, H., Yellapantula, S., de Frahan, M. T. H., Perry, B., 
Rood, J., Grout, R., & Day, M. (2021). Adaptive mesh based combustion simulations 
of direct fuel injection effects in a supersonic cavity flame-holder. 
Combustion and Flame, 232, 111531.

The runs from the paper are quite expensive. The baseline grid 
is 512 x 128 x 32 and simulations are performed with 3 levels (max_level=2).

First a non-reacting simulation is performed and then reacting cases for 
C1 and C2 are restarted from the non-reacting checkpoint file with 
inject_fuel=1 and appropriate injection location set by centx.

Further information is also available in the youtube video from 
APS gallery of fluid motion:
https://www.youtube.com/watch?v=QIaLfaethIw
