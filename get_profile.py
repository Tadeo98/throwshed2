import numpy as np
from throwshed2 import int_function

def get_profile(end_cell, SP, IH):
    dX = end_cell.GetX() - SP.GetX()
    dY = end_cell.GetY() - SP.GetY()
    try:
        Azimuth = np.arctan(dX / dY)
    except ZeroDivisionError:
        # for the case of dY being 0, making the division impossible
        if dX > 0:
            Azimuth = np.radians(90)
        else:
            Azimuth = np.radians(270)
            # azimuth needs to be recalculated accordingly to correct quadrant
    if dY > 0:
        if dX < 0:
            Azimuth += np.radians(360)
    elif dY < 0:
        Azimuth += np.radians(180)
    profile = [[0],[SP.GetZ()-IH]]
    s = 0
    cell_dist = ((SP.GetX() - end_cell.GetX()) ** 2 + (SP.GetY() - end_cell.GetY()) ** 2) ** (1 / 2)
    while True:
        s += 0.5
        X_compare_point = SP.GetX() + s * np.sin(Azimuth)
        Y_compare_point = SP.GetY() + s * np.cos(Azimuth)
        Z_compare_point = int_function(X_compare_point, Y_compare_point)
        profile[0].append(s)
        profile[1].append(Z_compare_point)
        if s > cell_dist:
            break
    return profile

"""
to use, I add following code after trajectory_set() function:
    import get_profile
    from matplotlib import pyplot as plt
    end_cell = ogr.Geometry(ogr.wkbPoint)
    end_cell.AddPoint(-488366, point_geom.GetY())
    profile = get_profile.get_profile(end_cell,SP,IH)
    plt.figure(figsize=(32, 18))
    plt.plot(profile[0], profile[1], '-', linewidth=1)
    mdti = TS.index((max(TS, key=lambda x: x[1][0][-1])))
    for i in range(len(TS)):
        plt.plot(TS[i][1][0], TS[i][1][1], '.-', linewidth=1)


    # print(TStemp)
    # plt.plot(TStemp[1][0], TStemp[1][1], '-', linewidth=1)
    plt.xlim(190, 230)
    plt.ylim(430, 485)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f'filename.png', dpi=300)
    plt.show()
    exit()
"""