#######################################################################
## THROWSHED ##
#######################################################################
# Numerical methods #
#######################################################################

"""
Module containing functions with numerical methods computing ballistic trajectories and returning lists of
coordinates. Euler() is for Euler differential method and Heun() is for Heun differential method.
"""

import os
import numpy as np
import atmospheres as atm
import drag_tables as dt

def interpolate_C_D(C_D_t, Mach):
    """
    Function that returns drag coefficient whether it is constant or it is interpolated from chosen drag to Mach function.
    C_D_t - drag coefficient type or table or constant, Mach - Mach number of the projectile
    """
    # drag to Mach function saved within module drag_tables or given by user (as set or as a constant)
    if type(C_D_t) == str:
        d2M = dt.dt_dir.get(C_D_t)
    else:
        d2M = C_D_t
    # interpolate the drag coefficient from list or set as constant
    if type(d2M) == list:
        for i in range(len(d2M[0])):
            if Mach < d2M[0][i]:
                C_D = d2M[1][i - 1] + (d2M[1][i] - d2M[1][i - 1]) / (d2M[0][i] - d2M[0][i - 1]) * (Mach - d2M[0][i - 1])
                break
    else:
        C_D = d2M
    return C_D

def euler2D(Y_0, AD, d, M, V, VA, T_0, atm_t, C_D_t, G, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH):
    """
    Euler differential method, generates trajectory 2D coordinates list.
    Metric units are used (m, °C, kg, °)
    Y_0 = initial altitude, AD - air density, d - diameter, M - mass, V - velocity, VA - vertical angle,
    T_0 - temperature at shooting site, atm_t - atmosphere type, C_D_t - drag coefficient type or table or constant,
    G - gravitational acceleration, TSS - trajectory segment size (width/length), TSD - trajectory segment dimension,
    IP - irregular projectile, CSA - cross-sectional area, PD - peak drag, PA - peak area,
    OD - oscillation distance, OF - oscillation frequency, DMINH - minimal height within DTM
    """
    # Cross-sectional area of the projectile (calculated from diameter if desired)
    if not IP:
        CSA = np.pi*d**2/4
    # list containing X,Y coordinates of trajectory points
    points = [[0.0], [Y_0]]
    # overall time
    T = 0
    # Curved Distance of Ballistic Trajectory
    CDBT = 0
    # incomplete drag formula - without drag coefficient and velocity (fixed auxiliary parameter)
    C_i = -AD * CSA / (2 * M)
    V_x, V_y = V * np.cos(VA), V * np.sin(VA)
    # set X element of velocity to 0 if the value is too small (for shooting angle -90 or 90 degrees when V_x is calculated as 1e-16)
    if V_x < 0.0000000001:
        V_x = 0
    # set trajectory segment width to given trajectory segment size
    TSW = TSS if V_x else 0
    # cycle going through all trajectory elements
    while True:
        # for the case of TSS being entered as length (not width), width has to be computed
        if TSD:
            TSW = TSS*np.cos(np.arctan(V_y/V_x)) if V_x else 0
        # X,Y steps added to points list (and Time increased)
        points[0].append(points[0][-1]+TSW)
        points[1].append(points[1][-1] + (V_y / V_x * TSW if V_x else TSS*V_y/abs(V_y)))
        # time and CDBT increases if there is a frequency
        if OF:
            T += TSW / V_x if V_x else TSS / abs(V_y)
            if TSD:
                CDBT += TSS
            else:
                CDBT += TSW/np.cos(np.arctan(V_y/V_x)) if V_x else TSS
        # when last height is less than minimal DEM height, cycle breaks and last values are reinterpolated into minimal DEM height (to prevent errors in extreme situations of further functions)
        if points[1][-1] <= DMINH:
            # if the shooting point is on the cell with minimal height of DEM, there will be only 2 points for trajectories starting with angle <= 0 and these points can't be the same, so this is the only exception where last points of trajectories are not recalculated (interpolated) to minimal DEM height
            if len(points[1]) == 2:
                pass
            # but normally the last segment passing the minimal DEM height is clipped by this height and the last point is interpolated to this height
            else:
                points[1][-1] = DMINH
                points[0][-1] = points[0][-2] + (points[1][-1] - points[1][-2]) * V_x / V_y
            break
        # velocity
        V = (V_x**2 + V_y**2)**(1/2)
        # temperature at relative height Y (above/below shooting site)
        Y = points[1][-1] - points[1][0]
        T_y = (T_0 + 273.15) * np.exp(-(atm.K[atm_t][0] + atm.K[atm_t][1]*Y)*Y) - 273.15
        # sound speed in air
        a_0 = atm.a_0_f[atm_t] * (T_y + 273.15) ** (1 / 2)
        # Mach number
        Mach = V/a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D = interpolate_C_D(C_D_t, Mach)
        # for oscillating/rotating objects (irregular projectiles, wobbling arrows) drag and area are adjusted according to frequency and max distance, after OD there is only regular drag and area if there's frequency at all
        if not OF:
            # no oscillation, element in [1/m] containing parameters with recalculated air density and the independent value is changed from time to range x
            C_D_x = C_i*C_D*V*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/(V_x if V_x else abs(V_y))
        else:
            if OD != 0 and CDBT >= OD:
                # no adjustments
                C_D_x = C_i * C_D * V * np.exp(-(atm.h[atm_t][0] + atm.h[atm_t][1] * Y) * Y) / (V_x if V_x else abs(V_y))
            else:
                # oscillation coefficient computed goniometrically from the current state of the projectile (ranges from regular position 0 to peak 1)
                OC = abs(np.sin((T*OF)%1*2*np.pi))
                # OC recalculated linearly according to ratio of the distance travelled to oscillation distance
                if OD:
                    OC *= (OD-CDBT)/OD
                # drag and area differences are adjusted accordingly
                DDIFF, ADIFF = (PD-C_D)*OC, (PA-CSA)*OC
                # C_D with added drag according to current state of the projectile
                C_D_c = C_D + DDIFF
                # C_i with added cross-sectional area according to current state of the projectile
                C_i_c = -AD * (CSA+ADIFF) / (2 * M)
                C_D_x = C_i_c*C_D_c*V*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/(V_x if V_x else abs(V_y))
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x = C_D_x*V_x
        F_y = C_D_x*V_y-G/(V_x if V_x else abs(V_y))
        # new velocity elements
        V_x, V_y = V_x + F_x*TSW, V_y + (F_y*TSW if V_x else F_y*TSS)
    return points

def euler3D(Y_0, AD, d, M, V, VA, W_x, W_z, T_0, Phi, Azi, atm_t, C_D_t, G, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH):
    """
    Euler differential method, generates trajectory 3D coordinates list.
    Metric units are used (m, °C, kg, °)
    Y_0 = initial altitude, AD - air density, d - diameter, M - mass, V - velocity, VA - vertical angle,
    W_x and W_z - wind velocity elements, T_0 - temperature at shooting site, Phi - latitude of firing site,
    Azi - azimut of fire clockwise from north, atm_t - atmosphere type,
    C_D_t - drag coefficient type or table or constant, G - gravitational acceleration,
    TSS - trajectory segment size (width/length), TSD - trajectory segment dimension, IP - irregular projectile,
    CSA - cross-sectional area, PD - peak drag, PA - peak area,
    OD - oscillation distance, OF - oscillation frequency, DMINH - minimal height within DTM
    """
    # Cross-sectional area of the projectile (calculated from diameter if desired)
    if not IP:
        CSA = np.pi*d**2/4
    Ome = 0.00007292
    # list containing X,Y,Z coordinates of trajectory points
    points = [[0.0], [Y_0], [0.0]]
    # overall time
    T = 0
    # Curved Distance of Ballistic Trajectory
    CDBT = 0
    # incomplete drag formula - without drag coefficient and velocity (fixed auxiliary parameter)
    C_i = -AD * CSA / (2 * M)
    V_x, V_y, V_z = V * np.cos(VA), V * np.sin(VA), 0
    # set trajectory segment width to given trajectory segment size
    TSW = TSS
    # cycle going through all trajectory elements
    while True:
        # for the case of TSS being entered as length (not width), width has to be computed
        if TSD:
            TSW = TSS*np.cos(np.arctan(V_y/(V_x**2+V_z**2)**(1/2)))
        # X,Y,Z steps added to points list
        points[0].append(points[0][-1]+TSW)
        points[1].append(points[1][-1] + V_y / V_x * TSW)
        points[2].append(points[2][-1] + V_z / V_x * TSW)
        # time and CDBT increases if there is a frequency
        if OF:
            T += TSW / V_x
            if TSD:
                CDBT += TSS
            else:
                CDBT += TSW / np.cos(np.arctan(V_y/(V_x**2+V_z**2)**(1/2)))
        # when last height is less than minimal DEM height, cycle breaks and last values are reinterpolated into minimal DEM height (to prevent errors in extreme situations of further functions)
        if points[1][-1] <= DMINH:
            # if the shooting point is on the cell with minimal height of DEM, there will be only 2 points for trajectories starting with angle <= 0 and these points can't be the same, so this is the only exception where last points of trajectories are not recalculated (interpolated) to minimal DEM height
            if len(points[1]) == 2:
                pass
            # but normally the last segment passing the minimal DEM height is clipped by this height and the last point is interpolated to this height
            else:
                points[1][-1] = DMINH
                points[0][-1] = points[0][-2] + (points[1][-1] - points[1][-2]) / V_y * V_x
                points[2][-1] = points[2][-2] + (points[1][-1] - points[1][-2]) / V_y * V_z
            break
        # velocity of projectile relative to wind speed
        V_r = ((V_x - W_x)**2 + V_y**2 + (V_z - W_z)**2)**(1/2)
        # temperature at relative height Y (above/below shooting site)
        Y = points[1][-1] - points[1][0]
        T_y = (T_0 + 273.15) * np.exp(-(atm.K[atm_t][0] + atm.K[atm_t][1]*Y)*Y) - 273.15
        # sound speed in air
        a_0 = atm.a_0_f[atm_t] * (T_y + 273.15) ** (1 / 2)
        # Mach number
        Mach = V_r/a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D = interpolate_C_D(C_D_t, Mach)
        # for oscillating/rotating objects (irregular projectiles, wobbling arrows) drag and area are adjusted according to frequency and max distance, after OD there is only regular drag and area if there's frequency at all
        if not OF:
            # no oscillation, element in [1/m] containing parameters with recalculated air density and the independent value is changed from time to range x
            C_D_x = C_i*C_D*V_r*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/V_x
        else:
            if OD != 0 and CDBT >= OD:
                # no adjustments
                C_D_x = C_i*C_D*V_r*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/V_x
            else:
                # oscillation coefficient computed goniometrically from the current state of the projectile (ranges from regular position 0 to peak 1)
                OC = abs(np.sin((T * OF) % 1 * 2 * np.pi))
                # OC recalculated linearly according to ratio of the distance travelled to oscillation distance
                if OD:
                    OC *= (OD - CDBT) / OD
                # drag and area differences are adjusted accordingly
                DDIFF, ADIFF = (PD - C_D) * OC, (PA - CSA) * OC
                # C_D with added drag according to current state of the projectile
                C_D_c = C_D + DDIFF
                # C_i with added cross-sectional area according to current state of the projectile
                C_i_c = -AD * (CSA + ADIFF) / (2 * M)
                C_D_x = C_i_c*C_D_c*V_r*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/V_x
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x = C_D_x*(V_x-W_x)+2*Ome*(-V_y*np.cos(Phi)*np.sin(Azi)-V_z*np.sin(Phi))/V_x
        F_y = C_D_x*V_y-G/V_x+2*Ome*(V_x*np.cos(Phi)*np.sin(Azi)+V_z*np.cos(Phi)*np.cos(Azi))/V_x
        F_z = C_D_x*(V_z-W_z)+2*Ome*(V_x*np.sin(Phi)-V_y*np.cos(Phi)*np.cos(Azi))/V_x
        # new velocity elements
        V_x, V_y, V_z = V_x + F_x*TSW, V_y + F_y*TSW, V_z + F_z*TSW
    return points

def heun2D(Y_0, AD, d, M, V, VA, T_0, atm_t, C_D_t, G, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH):
    """
    Heun differential method, generates trajectory 2D coordinates list.
    Metric units are used (m, °C, kg, ...)
    Y_0 = initial altitude, AD - air density, d - diameter, M - mass, V - velocity, VA - vertical angle,
    T_0 - temperature at shooting site, atm_t - atmosphere type, C_D_t - drag coefficient type or table or constant,
    G - gravitational acceleration, TSS - trajectory segment size (width/length), TSD - trajectory segment dimension,
    IP - irregular projectile, CSA - cross-sectional area, PD - peak drag, PA - peak area,
    OD - oscillation distance, OF - oscillation frequency, DMINH - minimal height within DTM
    """
    # just a note
    if OF:
        print("Heun method does not support oscillating projectiles.")
    # Cross-sectional area of the projectile (calculated from diameter if desired)
    if not IP:
        CSA = np.pi*d**2/4
    # list containing X,Y coordinates of trajectory points
    points = [[0.0], [Y_0]]
    # overall time
    T = 0
    # incomplete drag formula - without drag coefficient and velocity (fixed auxiliary parameter)
    C_i = -AD * CSA / (2 * M)
    V_x, V_y = V * np.cos(VA), V * np.sin(VA)
    # set X element of velocity to 0 if the value is too small (for shooting angle -90 or 90 degrees when V_x is calculated as 1e-16)
    if V_x < 0.0000000001:
        V_x = 0
    # set trajectory segment width to given trajectory segment size
    TSW = TSS if V_x else 0
    # set new trajectory segment width to given trajectory segment size
    TSW_n = TSS if V_x else 0
    # cycle going through all trajectory elements
    while True:
        # for the case of TSS being entered as length (not width), width has to be computed
        if TSD:
            TSW = TSS * np.cos(np.arctan(V_y / V_x)) if V_x else 0
        # velocity
        V = (V_x**2 + V_y**2)**(1/2)
        # temperature at relative height Y (above/below shooting site)
        Y = points[1][-1] - points[1][0]
        T_y = (T_0 + 273.15) * np.exp(-(atm.K[atm_t][0] + atm.K[atm_t][1]*Y)*Y) - 273.15
        # sound speed in air
        a_0 = atm.a_0_f[atm_t] * (T_y + 273.15) ** (1 / 2)
        # Mach number
        Mach = V/a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D = interpolate_C_D(C_D_t, Mach)
        # element in [1/m] containing parameters with recalculated air density and the independent value is changed from time to range x
        C_D_x = C_i*C_D*V*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/(V_x if V_x else abs(V_y))
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x = C_D_x*V_x
        F_y = C_D_x*V_y-G/(V_x if V_x else abs(V_y))
        # PREDICTOR - CORRECTOR APPLICATION, NEW PARAMETERS
        # new velocity elements
        V_x_n, V_y_n = V_x + F_x*TSW, V_y + (F_y*TSW if V_x else F_y*TSS)
        # new velocity
        V_n = (V_x_n ** 2 + V_y_n ** 2) ** (1 / 2)
        # Mach number
        Mach_n = V_n / a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D_n = interpolate_C_D(C_D_t, Mach_n)
        # drag coefficient with recalculated air density and the independent value is changed from time to range x
        C_D_x_n = C_i * C_D_n * V_n * np.exp(-(atm.h[atm_t][0] + atm.h[atm_t][1] * Y) * Y) / (V_x_n if V_x_n else abs(V_y_n))
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x_n = C_D_x_n * V_x_n
        F_y_n = C_D_x_n * V_y_n - G / (V_x_n if V_x_n else abs(V_y_n))
        # for the case of TSS being entered as length (not width), width for new (second) element has to be computed
        if TSD:
            TSW_n = TSS * np.cos(np.arctan(V_y_n / V_x_n)) if V_x_n else 0
        # recalculated new velocity elements with applied corrector
        V_x_c, V_y_c = V_x + (F_x*TSW+F_x_n*TSW_n)/2, V_y + ((F_y*TSW+F_y_n*TSW_n)/2 if V_x_n else F_y_n*TSS)
        # in point of trajectory, where the projectile starts to fall down, no corrector is applied, it is as in euler method (averaging values from these two elements creates undesired results)
        if V_y_c/V_y < 0:
            # X,Y,Z steps added to points list
            points[0].append(points[0][-1] + TSW)
            points[1].append(points[1][-1] + (V_y/V_x * TSW if V_x else TSS*V_y/abs(V_y)))
            # T += TSW/V_x
        # normally applied corrector
        else:
            # at last, TSW has to be recalculated according to the average of the velocities
            if TSD:
                TSW = TSS * np.cos(np.arctan((V_y + V_y_c) / (V_x + V_x_c))) if V_x and V_x_n else 0
            # X,Y,Z steps added to points list
            points[0].append(points[0][-1]+TSW)
            points[1].append(points[1][-1] + ((V_y + V_y_c) / (V_x + V_x_c) * TSW if V_x and V_x_n else TSS*(V_y + V_y_c)/abs((V_y + V_y_c))))
            #T += 2*TSW/(V_x + V_x_c)
        # when last height is less than minimal DEM height, cycle breaks and last values are reinterpolated into minimal DEM height (to prevent errors in extreme situations of further functions)
        if points[1][-1] <= DMINH:
            # if the shooting point is on the cell with minimal height of DEM, there will be only 2 points for trajectories starting with angle <= 0 and these points can't be the same, so this is the only exception where last points of trajectories are not recalculated (interpolated) to minimal DEM height
            if len(points[1]) == 2:
                pass
            # but normally the last segment passing the minimal DEM height is clipped by this height and the last point is interpolated to this height
            else:
                points[1][-1] = DMINH
                points[0][-1] = points[0][-2] + (points[1][-1] - points[1][-2]) * (V_x + V_x_c) / (V_y + V_y_c)
            break
        V_x, V_y = V_x_c, V_y_c
    return points

def heun3D(Y_0, AD, d, M, V, VA, W_x, W_z, T_0, Phi, Azi, atm_t, C_D_t, G, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH):
    """
    Heun differential method, generates trajectory 3D coordinates list.
    Metric units are used (m, °C, kg, ...)
    Y_0 = initial altitude, AD - air density, d - diameter, M - mass, V - velocity, VA - vertical angle,
    W_x and W_z - wind velocity elements, T_0 - temperature at shooting site, Phi - latitude of firing site,
    Azi - azimut of fire clockwise from north, atm_t - atmosphere type,
    C_D_t - drag coefficient type or table or constant, G - gravitational acceleration,
    TSS - trajectory segment size (width/length), TSD - trajectory segment dimension, IP - irregular projectile,
    CSA - cross-sectional area, PD - peak drag, PA - peak area,
    OD - oscillation distance, OF - oscillation frequency, DMINH - minimal height within DTM
    """
    # just a note
    if OF:
        print("Heun method does not support oscillating projectiles.")
    # Cross-sectional area of the projectile (calculated from diameter if desired)
    if not IP:
        CSA = np.pi*d**2/4
    # Angular velocity of Earth
    Ome = 0.00007292
    # list containing X,Y,Z coordinates of trajectory points
    points = [[0.0], [Y_0], [0.0]]
    # overall time
    T = 0
    # incomplete drag formula - without drag coefficient and velocity (fixed auxiliary parameter)
    C_i = -AD * CSA / (2 * M)
    V_x, V_y, V_z = V * np.cos(VA), V * np.sin(VA), 0
    # set trajectory segment width to given trajectory segment size
    TSW = TSS
    # set new trajectory segment width to given trajectory segment size
    TSW_n = TSS
    # cycle going through all trajectory elements
    while True:
        # for the case of TSS being entered as length (not width), width has to be computed
        if TSD:
            TSW = TSS * np.cos(np.arctan(V_y / (V_x ** 2 + V_z ** 2) ** (1 / 2)))
        # velocity of projectile relative to wind speed
        V_r = ((V_x - W_x)**2 + V_y**2 + (V_z - W_z)**2)**(1/2)
        # temperature at relative height Y (above/below shooting site)
        Y = points[1][-1] - points[1][0]
        T_y = (T_0 + 273.15) * np.exp(-(atm.K[atm_t][0] + atm.K[atm_t][1]*Y)*Y) - 273.15
        # sound speed in air
        a_0 = atm.a_0_f[atm_t] * (T_y + 273.15) ** (1 / 2)
        # Mach number
        Mach = V_r/a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D = interpolate_C_D(C_D_t, Mach)
        # element in [1/m] containing parameters with recalculated air density and the independent value is changed from time to range x
        C_D_x = C_i*C_D*V_r*np.exp(-(atm.h[atm_t][0]+atm.h[atm_t][1]*Y)*Y)/V_x
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x = C_D_x*(V_x-W_x)+2*Ome*(-V_y*np.cos(Phi)*np.sin(Azi)-V_z*np.sin(Phi))/V_x
        F_y = C_D_x*V_y-G/V_x+2*Ome*(V_x*np.cos(Phi)*np.sin(Azi)+V_z*np.cos(Phi)*np.cos(Azi))/V_x
        F_z = C_D_x*(V_z-W_z)+2*Ome*(V_x*np.sin(Phi)-V_y*np.cos(Phi)*np.cos(Azi))/V_x
        # PREDICTOR - CORRECTOR APPLICATION, NEW PARAMETERS
        # new velocity elements
        V_x_n, V_y_n, V_z_n = V_x + F_x*TSW, V_y + F_y*TSW, V_z + F_z*TSW
        # new relative velocity
        V_r_n = ((V_x_n - W_x) ** 2 + V_y_n ** 2 + (V_z_n - W_z) ** 2) ** (1 / 2)
        # Mach number
        Mach_n = V_r_n / a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D_n = interpolate_C_D(C_D_t, Mach_n)
        # drag coefficient with recalculated air density and the independent value is changed from time to range x
        C_D_x_n = C_i * C_D_n * V_r_n * np.exp(-(atm.h[atm_t][0] + atm.h[atm_t][1] * Y) * Y) / V_x_n
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x_n = C_D_x_n * (V_x_n - W_x)+2*Ome*(-V_y_n*np.cos(Phi)*np.sin(Azi)-V_z_n*np.sin(Phi))/V_x_n
        F_y_n = C_D_x_n * V_y_n - G / V_x_n+2*Ome*(V_x_n*np.cos(Phi)*np.sin(Azi)+V_z_n*np.cos(Phi)*np.cos(Azi))/V_x_n
        F_z_n = C_D_x_n * (V_z_n - W_z)+2*Ome*(V_x_n*np.sin(Phi)-V_y_n*np.cos(Phi)*np.cos(Azi))/V_x_n
        # for the case of TSS being entered as length (not width), width for new (second) element has to be computed
        if TSD:
            TSW_n = TSS * np.cos(np.arctan(V_y_n / (V_x_n ** 2 + V_z_n ** 2) ** (1 / 2)))
        # recalculated new velocity elements with applied corrector
        V_x_c, V_y_c, V_z_c = V_x + (F_x*TSW+F_x_n*TSW_n)/2, V_y + (F_y*TSW+F_y_n*TSW_n)/2, V_z + (F_z*TSW+F_z_n*TSW_n)/2
        # in point of trajectory, where the projectile starts to fall down, no corrector is applied, it is as in euler method (averaging values from these two elements creates undesired results)
        if V_y_c/V_y < 0:
            # X,Y,Z steps added to points list
            points[0].append(points[0][-1] + TSW)
            points[1].append(points[1][-1] + V_y/V_x * TSW)
            points[2].append(points[2][-1] + V_z/V_x * TSW)
            # T += TSW/V_x
        # normally applied corrector
        else:
            # at last, TSW has to be recalculated according to the average of the velocities
            if TSD:
                TSW = TSS * np.cos(np.arctan((V_y + V_y_c) / 2 / (((V_x + V_x_c) / 2) ** 2 + ((V_z + V_z_c) / 2) ** 2) ** (1 / 2)))
            # X,Y,Z steps added to points list
            points[0].append(points[0][-1] + TSW)
            points[1].append(points[1][-1] + (V_y + V_y_c) / (V_x + V_x_c) * TSW)
            points[2].append(points[2][-1] + (V_z + V_z_c) / (V_x + V_x_c) * TSW)
            # T += 2*TSW/(V_x + V_x_c)
        # when last height is less than minimal DEM height, cycle breaks and last values are reinterpolated into minimal DEM height (to prevent errors in extreme situations of further functions)
        if points[1][-1] <= DMINH:
            # if the shooting point is on the cell with minimal height of DEM, there will be only 2 points for trajectories starting with angle <= 0 and these points can't be the same, so this is the only exception where last points of trajectories are not recalculated (interpolated) to minimal DEM height
            if len(points[1]) == 2:
                pass
            # but normally the last segment passing the minimal DEM height is clipped by this height and the last point is interpolated to this height
            else:
                points[1][-1] = DMINH
                points[0][-1] = points[0][-2] + (points[1][-1] - points[1][-2]) / V_y * V_x
                points[2][-1] = points[2][-2] + (points[1][-1] - points[1][-2]) / V_y * V_z
            break
        V_x, V_y, V_z = V_x_c, V_y_c, V_z_c
    return points

def heun_i(Y_0, AD, d, M, V, VA, W_x, W_z, T_0, Phi, Azi, atm_t, C_D_t, G, TSW, DMINH):
    """
    Heun differential method, generates trajectory coordinates list.
    Imperial units are used (foot, °F, lb, inch...) as in original from McCoy's Modern Exterior Ballistics QBASIC script
    Y_0 = initial altitude, AD - air density, d - diameter, M - mass, V - velocity, VA - vertical angle,
    W_x and W_z - wind velocity elements, T_0 - temperature at shooting site, Phi - latitude of firing site,
    Azi - azimut of fire clockwise from north, atm_t - atmosphere type,
    C_D_t - drag coefficient type or table or constant, G - gravitational acceleration, TSW - trajectory segment width,
    DMINH - minimal height within DTM
    """
    # Cross-sectional area of the projectile
    CSA = np.pi * d ** 2 / 4
    # Angular velocity of Earth
    Ome = 0.00007292
    # list containing X,Y,Z coordinates of trajectory points
    points = [[0.0], [Y_0], [0.0]]
    # overall time
    T = 0
    # incomplete drag formula - without drag coefficient and velocity (fixed auxiliary parameter)
    C_i = -AD * CSA / (2 * M) / 144
    V_x, V_y, V_z = V * np.cos(VA), V * np.sin(VA), 0
    # cycle going through all trajectory elements
    while True:
        # velocity of projectile relative to wind speed
        V_r = ((V_x - W_x)**2 + V_y**2 + (V_z - W_z)**2)**(1/2)
        # temperature at relative height Y (above/below shooting site)
        Y = points[1][-1] - points[1][0]
        T_y = (T_0 + 459.67) * np.exp(-(atm.K_i[atm_t][0] + atm.K_i[atm_t][1]*Y)*Y) - 459.67
        # sound speed in air
        a_0 = atm.a_0_f_i[atm_t] * (T_y + 459.67)**(1/2)
        # Mach number
        Mach = V_r/a_0
        # interpolate Drag Coefficient (or can be constant)
        C_D = interpolate_C_D(C_D_t, Mach)
        # drag coefficient with recalculated air density and the independent value is changed from time to range x
        C_D_x = C_i*C_D*V_r*np.exp(-(atm.h_i[atm_t][0]+atm.h_i[atm_t][1]*Y)*Y)/V_x
        # elements in [1/s] simply called frequency, will be multiplied by the range step size
        F_x = C_D_x*(V_x-W_x)+2*Ome*(-V_y*np.cos(Phi)*np.sin(Azi)-V_z*np.sin(Phi))/V_x
        F_y = C_D_x*V_y-G/V_x+2*Ome*(V_x*np.cos(Phi)*np.sin(Azi)+V_z*np.cos(Phi)*np.cos(Azi))/V_x
        F_z = C_D_x*(V_z-W_z)+2*Ome*(V_x*np.sin(Phi)-V_y*np.cos(Phi)*np.cos(Azi))/V_x
        # PREDICTOR - CORRECTOR APPLICATION, NEW PARAMETERS
        # new velocity elements
        V_x_n, V_y_n, V_z_n = V_x + F_x*TSW, V_y + F_y*TSW, V_z + F_z*TSW
        # new relative velocity
        V_r_n = ((V_x_n - W_x) ** 2 + V_y_n ** 2 + (V_z_n - W_z) ** 2) ** (1 / 2)
        # cycle applying corrector within one element until velocity condition is met
        while True:
            # Mach number
            Mach_n = V_r_n / a_0
            # interpolate Drag Coefficient (or can be constant)
            C_D_n = interpolate_C_D(C_D_t, Mach_n)
            # drag coefficient with recalculated air density and the independent value is changed from time to range x
            C_D_x_n = C_i * C_D_n * V_r_n * np.exp(-(atm.h_i[atm_t][0] + atm.h_i[atm_t][1] * Y) * Y) / V_x_n
            # elements in [1/s] simply called frequency, will be multiplied by the range step size
            F_x_n = C_D_x_n * (V_x_n - W_x) + 2 * Ome * (-V_y_n * np.cos(Phi) * np.sin(Azi) - V_z_n * np.sin(Phi)) / V_x_n
            F_y_n = C_D_x_n * V_y_n - G / V_x_n + 2 * Ome * (V_x_n * np.cos(Phi) * np.sin(Azi) + V_z_n * np.cos(Phi) * np.cos(Azi)) / V_x_n
            F_z_n = C_D_x_n * (V_z_n - W_z) + 2 * Ome * (V_x_n * np.sin(Phi) - V_y_n * np.cos(Phi) * np.cos(Azi)) / V_x_n
            # recalculated new velocity elements with applied corrector
            V_x_c, V_y_c, V_z_c = V_x + (F_x+F_x_n)/2*TSW, V_y + (F_y+F_y_n)/2*TSW, V_z + (F_z+F_z_n)/2*TSW
            # relative velocity after applying corrector
            V_r_c = ((V_x_c - W_x) ** 2 + V_y_c ** 2 + (V_z_c - W_z) ** 2) ** (1 / 2)
            # velocity condition needs to be met, if not, new predictor-corrector will be applied, if yes corrector was applied successfully and the coordinates element's border are saved
            if abs((V_r_c-V_r_n)/V_r_c) > 0.00001:
                V_r_n = V_r_c
                V_x_n, V_y_n, V_z_n = V_x_c, V_y_c, V_z_c
            else:
                # X,Y,Z steps added to points list
                points[0].append(points[0][-1]+TSW)
                points[1].append(points[1][-1] + (V_y + V_y_c) / (V_x + V_x_c) * TSW)
                points[2].append(points[2][-1] + (V_z + V_z_c) / (V_x + V_x_c) * TSW)
                T += 2 * TSW / (V_x + V_x_c)
                break
        V_x, V_y, V_z = V_x_c, V_y_c, V_z_c
        if points[1][-1] <= DMINH:
            break
    return points