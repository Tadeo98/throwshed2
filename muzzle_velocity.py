#######################################################################
## MUZZLE VELOCITY ##
#######################################################################

"""
Allows to calculate muzzle velocity of cannons by formulae of internal velocity based on Benjamin Robins' work.
Parameters of cannon, cannonball and charge are required.
"""

import numpy as np

def main(formula, R, P_atm, p, m, eta, d, L):
    c = p*4/np.pi/d**2/eta #length of charge
    if formula == 1:
        v = (2*R*P_atm*p/m/eta*np.log(L/c))**(1/2)
    elif formula == 2:
        v = 587.7*(p/(m+p/3)*np.log(L/c))**(1/2)
    elif formula == 3:
        v = 606.9*(p/(m+p/3)*np.log(L/c))**(1/2)
    print('Muzzle velocity is:', np.round(v,2), 'm/s')

## SETTINGS
formula = 1 #type of the formula used, 1 = original Robins' model, 2 = refined formula for 18th century cannons, 3 = refined formula for 19th century cannons

## VARIABLES
R = 1600 #initial ratio of hot gas pressure to atmospheric pressure, 1500 for 18th century, 1600 for 19th century
P_atm = 101325 #air pressure [Pa]
p = 3.0801 #gunpowder charge weight [kg]
m = 8.1608 #cannonball weight [kg]
eta = 940 #density of gunpowder [kg/m^3]
d = 0.1364 #cannon barrel diameter [m]
L = 2.8376 #cannon barrel length [m]


if __name__ == "__main__":
    main(formula, R, P_atm, p, m, eta, d, L)