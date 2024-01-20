#######################################################################
## THROWSHED ##
#######################################################################
# Standard atmospheres #
#######################################################################

"""
Module containing variables and constants for Army Standard Metro (0) and ICAO (1) standard atmospheres to calculate
temperatures and air density at specific heights above shooting sites.
First field contains Army Standard Metro values and second field contains ICAO values.
"""

# CONSTANTS [Army Standard Metro, ICAO]

# Temperature-altitude decay factor
K = [[0.000019734, 0.0], [0.000022500, 0.000000000298806]]
# Speed of sound in air factor
a_0_f = [20.1186, 20.05]
# Air density-altitude decay factor
h = [[0.00010361, 0.0], [0.00009600, 0.00000000108]]

K_i = [[0.000006015, 0.0], [0.000006858, 0.00000000002776]]
# Speed of sound in air factor
a_0_f_i = [49.19, 49.0223]
# Air density-altitude decay factor
h_i = [[0.00003158, 0.0], [0.00002926, 0.0000000001]]

#3.28084