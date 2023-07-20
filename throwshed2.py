#######################################################################
## THROWSHED ##
#######################################################################
import os
import numpy as np
from osgeo import gdal, ogr, osr


#######################################################################
## FUNCTIONS

def main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, use_viewshed, EPSG,
         cumulative_throwshed, initial_height, initial_velocity, drag_coefficient, cross_sectional_area, mass,
         eyes_height, target_height, wall_height, constant, area_addition, wobble_distance, band_number=1,
         interpolation=1, alpha_min=-90.0, alpha_max=90.0, gravitational_acceleration=-9.81, air_density=1.225, dalpha=5,
         trajectory_segment_width=None):
    """Just main function with controls, global variables settings and triggers to other functions"""
    # making sure the vertical angle has correct range
    if alpha_max < alpha_min:
        print("Minimal vertical shooting angle higher than maximal.")
        exit()
    # Global variables
    global SRS, DP, PLP, LLP, TOF, TF, UV, CV, BN, INT, BF, VF, TSW, AL, DDS, DB, DA, DGT, DMINH, DMAXH, IH, IV, DC, \
        CSA, M, CONST, AA, WD, GA, AD, EH, TH, NDV, TA
    # CRS and other variable definition
    SRS = osr.SpatialReference()
    SRS.ImportFromEPSG(EPSG)
    TOF, TF, UV, CV, BN, INT, VF, IH, IV, DC, CSA, M, CONST, AA, WD, GA, AD, EH, TH = throwshed_output_folder,\
        throwshed_file, use_viewshed, cumulative_throwshed, band_number, interpolation, 'viewshed', initial_height, \
        initial_velocity, drag_coefficient, cross_sectional_area, mass, constant, area_addition, wobble_distance, \
        gravitational_acceleration, air_density, eyes_height, target_height
    # get DEM data and assign them as global (and referencing datasource)
    DDS, DB, DA, DGT, NDV = get_raster_from_file(dem_path)
    # assign trajectory segment width
    TSW = np.min(np.abs([DGT[1],DGT[5]])) if trajectory_segment_width == None else trajectory_segment_width
    # obtain list of point geometries (and all referencing data it's dependent on)
    point_layer_ds, point_layer, point_feature_list, point_geom_list = get_geom_list_from_file(point_layer_path)
    # burn lines as obstacles into DEM, creating new DEM (Digital terrain model -> Digital surface model)
    if line_layer_path:
        burn_obstacles(line_layer_path, wall_height)
    # get minimum and maximum DEM height
    DMINH, DMAXH = get_min_max_height()
    # obtain list of vertical angles
    AL = np.arange(np.radians(alpha_min), np.radians(alpha_max + dalpha), np.radians(dalpha))
    AL[-1] = np.radians(alpha_max) # in case of last angle being larger than 90°
    # throwshed array containing zeroes at first, will be edited later, dimensions same as dimensions of DEM raster
    TA = [np.zeros((DA.shape[0], DA.shape[1]), np.int16), np.zeros((DA.shape[0], DA.shape[1]), np.int16)]
    # cycle calculating throwshed for each point
    for i, point_geom in enumerate(point_geom_list):
        # compute throwshed for 1 point
        throwshed(point_geom, i)
        #plot_trajectory()
        break
    # finally, array is written into band of output throwshed raster
    create_raster_file(TF, TA, 1, gdal.GDT_Int16)
    DDS = DB = DA = DGT = NDV = None



def get_raster_from_file(file_path):
    """Get DEM datasource, band, array and geotransformation data and nodata value"""
    # import DEM datasource
    dem_ds = gdal.Open(file_path)
    # select band
    dem_band = dem_ds.GetRasterBand(BN)
    # DEM cell values into array
    dem_array = dem_band.ReadAsArray()
    # transformation data describing DEM
    dem_gt = dem_ds.GetGeoTransform()
    # nodata value
    no_data_value = dem_band.GetNoDataValue()
    return dem_ds, dem_band, dem_array, dem_gt, no_data_value

def get_min_max_height():
    """Get minimum and maximum DEM height"""
    # sometimes the function returns None, therefore statistics need to be calculated first
    if DB.GetMinimum() == None or DB.GetMaximum() == None:
        DB.ComputeStatistics(0)
    return DB.GetMinimum(), DB.GetMaximum()

def get_geom_list_from_file(file_path):
    """Get features' geometries from vector layer file"""
    # import vector layer datasource
    layer_ds = ogr.Open(file_path, 0)  # 1 = editing, 0 = read only. Datasource
    # vector layer
    layer = layer_ds.GetLayer()
    # list of features, then geometries
    feature_list = [layer.GetFeature(number) for number in range(0, layer.GetFeatureCount())]
    geom_list = [feature.GetGeometryRef() for feature in feature_list]
    return layer_ds, layer, feature_list, geom_list

def burn_obstacles(line_layer_path, wall_height):
    """Lines that resemble walls/obstacles are burnt into DEM"""
    # get list of line geometries
    line_geom_list = get_geom_list_from_file(line_layer_path)
    # calculate minimal buffer distance, at which the obstacle will always be respected
    buffer_dist = TSW / 2 + (DGT[1] ** 2 + DGT[5] ** 2) ** (1 / 2) / 2 - (TSW / 2) / (
                DGT[1] ** 2 + DGT[5] ** 2) ** (1 / 2) % 1
    # 1st buffer is created and in case of multiple lines the cycle is utilized to create and unite buffers
    buffer_geom = line_geom_list[0].Buffer(buffer_dist, 1)
    if len(line_geom_list) > 1:
        for line_geom in line_geom_list[1:]:
            buffer_geom = buffer_geom.Union(line_geom.Buffer(buffer_dist, 1))
    # buffer file temporary name
    BFT = "buffer_temp"
    # create layer for buffer
    buffer_outlayer = create_outlayer(BFT, buffer_geom)
    # create buffer datasource
    buffer_ds, buffer_band = create_raster_file(BFT, None, 0, gdal.GDT_Float32)
    # Buffer polygon is rasterized
    gdal.RasterizeLayer(buffer_ds, [1], buffer_outlayer, burn_values=[wall_height])
    # Sum of initial dem and buffer rasters
    buffer_array = buffer_band.ReadAsArray()
    global DA, DB
    DA = np.add(DA, buffer_array)
    buffer_ds = buffer_outlayer = None
    # create new raster datasource for DEM (DSM)
    dem_ds, DB = create_raster_file(BFT, [DA], 1, gdal.GDT_Float32)
    dem_ds = None
    # delete all temporary files
    remove_temp_files(BFT)


def create_outlayer(layer_name, geom):
    """Creates file for output layer, which will contain new features"""
    # create driver and output data source
    outds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(TOF + "\\" + layer_name + ".shp")
    # create output layer
    outlayer = outds.CreateLayer(layer_name, SRS)
    # feature definition and setting
    feature = ogr.Feature(outlayer.GetLayerDefn())
    feature.SetGeometry(geom)
    # assign feature into output layer
    outlayer.CreateFeature(feature)
    outds = None
    return outlayer

def create_raster_file(raster_name, dem_array_list, method, GDT):
    """Creates raster file. Method 0 returns empty datasource and band. Method 1 returns datasource and band with written array"""
    # create driver and output data source
    outds = gdal.GetDriverByName('GTiff').Create(TOF + "\\" + raster_name + ".tif", xsize=DA.shape[1],
                                                 ysize=DA.shape[0], bands=len(dem_array_list), eType=GDT)
    # assign geotransformation, projection, band and nodata settings
    outds.SetGeoTransform(DGT)
    outds.SetProjection(SRS.ExportToWkt())
    for i, dem_array in enumerate(dem_array_list):
        raster_band = outds.GetRasterBand(i+1)
        raster_band.SetNoDataValue(NDV)
        if method:
            raster_band.WriteArray(dem_array)
    return outds, raster_band

def remove_temp_files(temp_file):
    """Deletes all temporary files with assigned name"""
    for format in [file.split('.')[1] for file in os.listdir(TOF) if file.split('.')[0] == temp_file]:
        os.remove(TOF + '\\' + temp_file + "." + format)

def throwshed(point_geom, k):
    """Calculates throwshed for 1 point"""
    global SP
    # create shooting point with Z coordinate that is interpolated from DEM
    SP = ogr.Geometry(ogr.wkbPoint)
    SP.AddPoint(point_geom.GetX(), point_geom.GetY(), float(int_function(point_geom.GetX(), point_geom.GetY()))+IH)
    # generate set of trajectories for vertical angle range with basic step
    trajectory_simple_set()
    # insert trajectories between those from simple set, to ensure throwshed's edge accuracy
    trajectory_initial_set()

    plot_trajectory()
    exit()

    # define Ascending and Descending Trajectory Fields
    create_trajectory_fields()
    # Assign values into arrays of 2 bands (ATF and DTF)
    assign_values_to_throwshed(k)

def int_function(X, Y):
    """Interpolates height of point from DEM cells"""
    # nearest neighbour
    if INT == 0:
        column = round(np.abs((X - (DGT[0] + DGT[1] / 2)) / DGT[1]))
        row = round(np.abs((Y - (DGT[3] + DGT[5] / 2)) / DGT[5]))
        return DA[row][column]
    # bilinear
    else:
        left_column = int(np.abs((X - (DGT[0] + DGT[1] / 2)) / DGT[1])) # index of left column in set of four cells
        upper_row = int(np.abs((Y - (DGT[3] + DGT[5] / 2)) / DGT[5])) # index of the upper row in set of four cells
        X_left_cell = DGT[0] + DGT[1] / 2 + left_column * DGT[1] # X coordinate of left cells
        Y_lower_cell = DGT[3] + DGT[5] / 2 + (upper_row + 1) * DGT[5] # Y coordinate of lower cells
        H_1 = DA[upper_row][left_column]  # height of upper left cell
        H_2 = DA[upper_row][left_column + 1]  # height of upper right cell
        H_3 = DA[upper_row + 1][left_column]  # height of lower left cell
        H_4 = DA[upper_row + 1][left_column + 1]  # height of lower right cell
        H_int_1 = ((X - X_left_cell) * (H_4 - H_3)) / (np.abs(DGT[1])) + H_3  # interpolated height among lower cells
        H_int_2 = ((X - X_left_cell) * (H_2 - H_1)) / (np.abs(DGT[1])) + H_1  # interpolated height among upper cells
        return ((Y - Y_lower_cell) * (H_int_2 - H_int_1)) / (np.abs(DGT[5])) + H_int_1

def trajectory_simple_set():
    """Generate set of trajectories for vertical angle range with basic step"""
    # trajectory dictionary, that will contain all trajectories, their initial shooting angle etc.
    global TS
    # element in list begins with alpha value and continues with list of x and y coords lists in one trajectory
    TS = [[alpha, generate_trajectory(alpha)] for alpha in AL]


def generate_trajectory(alpha):
    """Generates trajectory from input parameters"""
    # initial drag
    d = -AD * IV ** 2 * DC * CONST * (CSA + AA) / (2 * M)
    # initial projectile velocity
    V = IV
    # list of all trajectory points' coordinates
    points = [[0.0], [SP.GetZ()]]
    # velocity x and y elements
    V_x = V * np.cos(alpha)
    V_y = V * np.sin(alpha)
    # drag x and y elements
    d_x = d * np.cos(alpha)
    d_y = d * np.sin(alpha)
    # time step is set to value approximate to half of cell size
    dt = TSW / V / 2
    # x and y steps
    dX = V_x * dt + d_x / 2 * dt ** 2
    dY = V_y * dt + (d_y + GA) / 2 * dt ** 2
    # cycle calculating every new trajectory point one by one
    while True:
        # coords
        points[0].append(points[0][-1] + dX)
        points[1].append(points[1][-1] + dY)
        # when last height is less than minimal DEM height, cycle breaks and last values are reinterpolated into minimal DEM height (to prevent errors in extreme situations of further functions)
        if points[1][-1] <= DMINH:
            # if the shooting point is on the cell with minimal height of DEM, there will be only 2 points for trajectories starting with angle <= 0 and these points can't be the same, so this is the only exception where last points of trajectories are not recalculated (interpolated) to minimal DEM height
            if len(points[1]) == 2:
                pass
            # but normally the last segment passing the minimal DEM height is clipped by this height and the last point is interpolated to this height
            else:
                points[1][-1] = DMINH
                points[0][-1] = points[0][-2] + (points[1][-1] - points[1][-2])/np.tan(alpha)
            break
        # new vertical angle
        alpha = np.arctan(dY / dX)
        # new velocity
        V = ((dX / dt) ** 2 + (dY / dt) ** 2) ** (1 / 2)
        # new drag
        if (points[0][-1] ** 2 + (points[1][-1] - SP.GetZ()) ** 2) ** (1 / 2) < WD:
            d = -AD * V ** 2 * DC * CONST * (CSA + AA) / (2 * M)
        else:
            d = -AD * V ** 2 * DC * CSA / (2 * M)
        # new drag x and y elements
        d_x = d * np.cos(alpha)
        d_y = d * np.sin(alpha)
        # new velocity x and y elements
        V_x += d_x * dt
        V_y += (d_y + GA) * dt
        # time step is recalculated to value approximate to half of cell size
        dt = TSW / V / 2
        # new x and y steps
        dX = V_x * dt + d_x / 2 * dt ** 2
        dY = V_y * dt + (d_y + GA) / 2 * dt ** 2
    return points

def trajectory_initial_set():
    """Calculates and inserts trajectories between those in simple set, to make it denser and to ensure throwshed's
    edge accuracy. Calculates and returns trajectory envelope points list (its useful section)."""
    global TS, envelope
    # new and previous trajectory end x, first ones are random, just to make sure the cycle does not stop immediately
    ntex = [(max(TS, key=lambda x: x[1][0][-1])[1][0][-1]+TSW)*2,(max(TS, key=lambda x: x[1][0][-1])[1][0][-1]+TSW)*3]
    # cycle that finds the furthest possible trajectory for minimal DEM height respecting the edge accuracy
    while round(np.abs(ntex[0] - ntex[1])/TSW):
        #most distant trajectory index
        mdti = TS.index((max(TS, key=lambda x: x[1][0][-1])))
        # adds new trajectories before and after current furthest trajectory
        if mdti != 0 and mdti != len(TS)-1:
            for new_alpha in [(TS[mdti+1][0] - TS[mdti][0]) / 2 + TS[mdti][0], (TS[mdti][0] - TS[mdti-1][0]) / 2 + TS[mdti-1][0]]:
                TS.append([new_alpha, generate_trajectory(new_alpha)])
        # for furthest trajectory that is also the first or last one, only one trajectory is added accordingly
        elif mdti == len(TS)-1:
            new_alpha = (TS[mdti][0] - TS[mdti-1][0]) / 2 + TS[mdti-1][0]
            TS.append([new_alpha, generate_trajectory(new_alpha)])
            mdti += 1 #this is just so that mdti gets higher like length of TS does (to cope with length of TS getting bigger, these 2 are compared in ntex)
        elif mdti == 0:
            new_alpha = (TS[mdti+1][0] - TS[mdti][0]) / 2 + TS[mdti][0]
            TS.append([new_alpha, generate_trajectory(new_alpha)])
        ntex = [max(TS[-1][1][0][-1],TS[-2][1][0][-1]) if mdti != 0 and mdti != len(TS)-1 else TS[-1][1][0][-1], ntex[0]]
        TS.sort(key=lambda x: x[0])

    # initial trajectory index (starting will be the trajectory with furthest reach)
    iti = TS.index((max(TS, key=lambda x: x[1][0][-1])))
    # function ends after finding out the last trajectory is the one with furthest reach (rest of the code is not applicable) and empty envelope is returned, otherwise trajectory set will get denser with following code
    if iti == len(TS)-1:
        envelope = [[], []]
        return
    # envelope needs starting point before adding more points to it
    envelope = [[TS[iti][1][0][-1]], [TS[iti][1][1][-1]]]
    # X and Y Inner Interection from Previous Cycle will be set to furthest point during first cycle
    XIIPR, YIIPR = envelope[0][0], envelope[1][0]
    # cycle that inserts trajectories between furthest trajectory at minimal DEM height and maximal DEM height
    try:
        while True:
            # generate new trajectory in between actual and following
            new_alpha = abs(TS[iti][0] - TS[iti + 1][0]) / 2 + TS[iti][0]
            TS.insert(iti+1, [new_alpha, generate_trajectory(new_alpha)])
            # in case of added trajectory having further reach than previously furthest trajectory (actual)
            if TS.index((max(TS, key=lambda x: x[1][0][-1]))) == iti + 1:
                iti += 1
                # even starting point of envelope needs to be updated
                envelope = [[TS[iti][1][0][-1]], [TS[iti][1][1][-1]]]
                XIIPR, YIIPR = envelope[0][0], envelope[1][0]
                continue

            # get intersection of actual and following (newly created) trajectory (Right Outer Intersection)
            XROI, YROI = intersection_of_trajectories(iti, iti+1)
            # get intersection of following and 2. following trajectory (Left Outer Intersection)
            XLOI, YLOI = intersection_of_trajectories(iti+1, iti+2)
            # get intersection of actual and 2. following trajectory (Inner Intersection), first one is calculated, the rest will be reused from outer intersections
            XII, YII = intersection_of_trajectories(iti, iti+2)

            # if the last trajectory incorporated in cycle is the last one from the net with shooting angle value of 90 degrees.
            # This is because with 90 degrees trajectory left outer and inner intersections will be the same, which would lead to undesired behaviour
            if TS[iti+2][0] == np.radians(90):
                # if X coordinate of highest point in newly created trajectory is less than TSW, last possible area of the net was made dense enough
                if not np.floor(TS[iti+1][1][0][TS[iti+1][1][1].index((max(TS[iti+1][1][1])))]/TSW):
                    # Initial Trajectory Reversed list (X and Y coords)
                    ITR = [TS[iti][1][0][-1::-1], TS[iti][1][1][-1::-1]]
                    # update envelope with points from initial trajectory of last cycle
                    update_envelope(ITR, XIIPR, XROI, YROI)
                    # Last Inserted Trajectory Reversed list (X and Y coords)
                    LITR = [TS[iti+1][1][0][-1::-1], TS[iti+1][1][1][-1::-1]]
                    # last points from last inserted trajectory and highest point of last trajectory are appended to envelope
                    for x, y in zip(LITR[0], LITR[1]):
                        if y > YROI:
                            envelope[0].append(x)
                            envelope[1].append(y)
                        if x == TS[iti+1][1][0][TS[iti+1][1][1].index((max(TS[iti+1][1][1])))]:
                            break
                    envelope[0].append(0)
                    envelope[1].append(max(TS[-1][1][1]))
                    break
                # if not dense enough, density will be accomplished with new iteration
                else:
                    continue

            # coordinates of midpoint on line between outer intersections (Outer Midpoint)
            XOM, YOM = (XROI + XLOI) / 2, (YROI + YLOI) / 2
            # following trajectory reversed list (X and Y coords)
            FTR = [TS[iti + 1][1][0][-1::-1], TS[iti + 1][1][1][-1::-1]]
            # finds intersection of arc distance and particular segment on the arc
            for i in range(len(FTR[0])-1):
                XA, YA = calculate_intersection(XOM, YOM, XII, YII, FTR[0][i], FTR[1][i], FTR[0][i+1], FTR[1][i+1])
                # checking if the intersection is really on the segment, if so, the Distance of Arc from Inner Intersection is calculated
                if FTR[0][i] >= XA >= FTR[0][i+1]:
                    DAII = ((XA-XII)**2 + (YA-YII)**2)**(1/2)
                    break
            # controls - compare horizontal distance of intersections and distance of arc from inner intersection
            if round(np.abs(XROI - XLOI)/TSW) and round(DAII/TSW):
                continue
            else:
                # with each shooting point the amount of these inserted auxiliary trajectories would almost double which could create pointless amount of trajectories
                del TS[iti+1]
                # initial trajectory reversed list (X and Y coords)
                ITR = [TS[iti][1][0][-1::-1], TS[iti][1][1][-1::-1]]
                # update envelope
                update_envelope(ITR, XIIPR, XII, YII)
                # previous intersection for next cycle is assigned
                XIIPR, YIIPR = XII, YII
                # at least one of the conditions was met and the cycle can jump to next initial trajectory
                iti += 1
                # if the cycle comes to last trajectory, it breaks as there is no following trajectory
                if TS[iti][0] == AL[-1]:
                    break
    except Exception as e:
        print(e)

def intersection_of_trajectories(t1i,t2i):
    """Looks for intersection between two trajectories and returns its X and Y coordinates.
    t1i and t2i are indexes of first and second trajectory between which the intersection is sought."""
    # reversed lists of trajectories' coordinates as the algorithm starts from end points, TX1 = X coordinates of 1. trajectory
    T1X, T1Y = TS[t1i][1][0][-1::-1], TS[t1i][1][1][-1::-1]
    T2X, T2Y = TS[t2i][1][0][-1::-1], TS[t2i][1][1][-1::-1]
    # X and Y coords of intersection (to be compared e.g. with XPI), i2s stands for radius around i2 (or index)
    XI = YI = i2s = False
    # following trajectory segment radius where the intersection will be sought
    ftsr = [0, 1, -1, 2, -2]
    # cycle that starts comparing coords of actual trajectory points from the end
    for i1 in range(1, len(T1X)):
        # cycle that starts comparing coords of following trajectory points from the end
        for i2 in range(1, len(T2X)):
            if T2Y[i2] > T1Y[i1] and T1X[i1] < T2X[i2 - 1]:
                # when potentially intersecting segment of following trajectory is found, because of rare situations its 2 following and preceding segments have to be assessed
                for i2s in ftsr:
                    XI, YI = calculate_intersection(T1X[i1 - 1], T1Y[i1 - 1], T1X[i1], T1Y[i1],
                                               T2X[i2 - 1 + i2s], T2Y[i2 - 1 + i2s], T2X[i2 + i2s],
                                               T2Y[i2 + i2s])
                    # making sure the intersection is between existing segments, not on their extension
                    if T1X[i1 - 1] >= XI >= T1X[i1] and T2X[i2 - 1 + i2s] >= XI >= T2X[i2 + i2s]:
                        break
                    XI = YI = False
            if XI or i2s == ftsr[-1]:
                # i2s makes sure that if the intersection is not found in the 2 segment radius of following trajectory, index for segment of actual trajectory has to increase
                i2s = False
                break
        if XI:
            break
    return XI, YI

def calculate_intersection(x1, y1, x2, y2, x3, y3, x4, y4):
    """Calculates intersection of 2 line segments and returns its X and Y coordinates"""
    x1, y1, x2, y2, x3, y3, x4, y4 = x1, y1, x2, y2, x3, y3, x4, y4
    XI = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) *
        (x3 - x4))
    YI = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) *
        (x3 - x4))
    return XI, YI

def update_envelope(ITR, XIIPR, XII, YII):
    """Updates envelope with parts of trajectories."""
    global envelope
    # on reversed trajectory, index of first point that will be incorporated within the envelope is found
    for i in range(len(ITR[0]) - 1):
        if ITR[0][i] >= XIIPR >= ITR[0][i + 1]:
            # Envelope Part Starting Point Index
            EPSPI = i + 1
            break
    # on reversed trajectory, index of last point that will be incorporated within the envelope is found
    for i in range(EPSPI - 1, len(ITR[0]) - 1):
        if ITR[0][i] >= XII >= ITR[0][i + 1]:
            # Envelope Part Ending Point Index
            EPEPI = i
            break
    # condition for rare situation where both inner intersections could fall within one segment of trajectory
    if EPEPI >= EPSPI:
        # envelope is updated with all points between starting and ending point
        for i in range(EPSPI, EPEPI + 1):
            envelope[0].append(ITR[0][i])
            envelope[1].append(ITR[1][i])
    # lastly, envelope is updated with the inner intersection point
    envelope[0].append(XII)
    envelope[1].append(YII)

def create_trajectory_fields():
    """Creates ATF - Ascending Trajectory Field and DTF - Descending Trajectory Field. Lists are made into polygons."""
    global ATF_polygon, DTF_polygon
    ATF_polygon, DTF_polygon = None, None
    ATF, DTF = [[], []], [[], []]
    # in this case DTF does not exist
    if len(envelope[0]) == 0:
        ATF[0] = TS[0][1][0] + TS[-1][1][0][-1::-1]
        ATF[1] = TS[0][1][1] + TS[-1][1][1][-1::-1]
        # create polygon out of the list
        ATF_polygon = create_polygon_from_coords_list(ATF)
    # in this case both ATF and DTF are created
    else:
        # find out index of coords from last trajectory, where envelope connects with it
        for i in range(len(TS[-1][1][0])-1):
            # for last trajectory with shooting angle 90 degrees this stops immediately, to ATF is added first point of last trajectory, which still creates correct polygon
            if TS[-1][1][0][i] <= envelope[0][-1] <= TS[-1][1][0][i+1]:
                break
        ATF[0] = TS[0][1][0] + envelope[0] + TS[-1][1][0][i::-1]
        ATF[1] = TS[0][1][1] + envelope[1] + TS[-1][1][1][i::-1]
        # but for DTF with last trajectory shooting angle being 90 degrees, i is edited so that only last point down at min DEM height is added to polygon
        i = -2 if TS[-1][0] == AL[-1] else i
        DTF[0] = envelope[0] + TS[-1][1][0][i+1:] + envelope[0][:1]
        DTF[1] = envelope[1] + TS[-1][1][1][i+1:] + envelope[1][:1]
        # create polygons out of the lists
        ATF_polygon = create_polygon_from_coords_list(ATF)
        DTF_polygon = create_polygon_from_coords_list(DTF)

def create_polygon_from_coords_list(x_y_list):
    """Creates ring from list of X and Y coordinates, then uses ring to create polygon which is returned."""
    # create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in zip(x_y_list[0], x_y_list[1]):
        ring.AddPoint(x, y)
    # create polygon
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    return polygon

def assign_values_to_throwshed(k):
    """Assigns/adds values into throwshed arrays."""
    # cycle going through every single cell of DEM
    for i in range(DA.shape[0]):
        print(i)
        for j in range(DA.shape[1]):
            # with multiple shooting points nodata value can already be assigned to the cell, therefore the algorithm jumps to following cell
            if TA[0][i][j] == NDV:
                continue
            # nodata value is assigned to both arrays for both bands (ATF and DTF)
            if k == 0 and DA[i][j] == NDV:
                TA[0][i][j] = NDV
                TA[1][i][j] = NDV
                continue
            # calculate coordinates of cell's middle point and its horizontal distance from shooting point
            X_coor_cell = DGT[0]+(j+1/2)*DGT[1]
            Y_coor_cell = DGT[3]+(i+1/2)*DGT[5]
            cell_distance = ((SP.GetX() - X_coor_cell)**2 + (SP.GetY()-Y_coor_cell)**2)**(1/2)
            # create cell point with relative coordinates in the plane of trajectories and find out whether it's within the field, if so, further actions are conducted
            relative_cell = ogr.Geometry(ogr.wkbPoint)
            relative_cell.AddPoint(cell_distance, float(DA[i][j]))
            # call function to find cell intersecting trajectory and to determine whether the cell is reachable without any obstacles
            if relative_cell.Within(ATF_polygon):
                if find_intersecting_trajectory(1, -1, -1, -1, relative_cell):
                    TA[0][i][j] += 1
            # can be None
            if DTF_polygon:
                if relative_cell.Within(DTF_polygon):
                    if find_intersecting_trajectory(-1, len(TS), len(TS), len(TS), relative_cell):
                        TA[1][i][j] += 1

def find_intersecting_trajectory(step, i, i1, i2, relative_cell):
    """Finds trajectory that intersects the cell (or is close enough, within allowed distance). For ATF cycle
    increments from start to end of trajectory set and viceversa for DTF. Returns True if the cell is accessible
    or False if the cell is not accessible without any obstacles - this is determined further function."""
    # at first, 2 surrounding trajectories are found by making polygon out of them and asking whether the cell lies within
    while True:
        i += step
        # create polygon specially for the last two trajectories when looking for ascending or first two when looking for descending trajectory if the last one has shooting angle 90 degrees
        if i == len(TS) - 2 and step == 1 and TS[-1][0] == np.radians(90):
            polygon = create_polygon_from_coords_list([TS[i][1][0][:TS[i][1][1].index((max(TS[i][1][1])))+1] + TS[i+step][1][0][TS[i+step][1][1].index((max(TS[i+step][1][1])))::-1], TS[i][1][1][:TS[i][1][1].index((max(TS[i][1][1])))+1] + TS[i + step][1][1][TS[i+step][1][1].index((max(TS[i+step][1][1])))::-1]])
        elif i == len(TS) - 1 and step == -1 and TS[-1][0] == np.radians(90):
            polygon = create_polygon_from_coords_list([TS[i][1][0][-1:TS[i][1][1].index((max(TS[i][1][1])))+1:-1] + TS[i+step][1][0][TS[i+step][1][1].index((max(TS[i+step][1][1]))):] + TS[i][1][0][-1:], TS[i][1][1][-1:TS[i][1][1].index((max(TS[i][1][1])))+1:-1] + TS[i + step][1][1][TS[i+step][1][1].index((max(TS[i+step][1][1]))):] + TS[i][1][1][-1:]])
        # or create polygon regularly
        else:
            polygon = create_polygon_from_coords_list([TS[i][1][0] + TS[i + step][1][0][-1::-1], TS[i][1][1] + TS[i + step][1][1][-1::-1]])
        # if the cell lies within, normals are computed to assess the smallest perpendicular distance from trajectory to cell
        if relative_cell.Within(polygon):
            # this is just to recycle normal1 computed in previous cycle that will be the same if first half of previous polygon is the one the cell lies within
            if i1 - i:
                # compute normal from previous trajectory segment to cell
                normal1, j = compute_normal(i, relative_cell.GetX(), relative_cell.GetY())
            # auxiliary index i1 to find out if the normal1 was already computed, will be used in next iteration, works with both directions of incrementing
            if step == -1:
                i1 = i + 1
            else:
                i1 = i
            # if cell is not close enough to trajectory to consider it as piercing trajectory, following trajectory is tested. If it is close enough, cycle breaks with i being the index of intersecting trajectory
            if np.round(normal1 / TSW):
                # this is just to recycle normal2 computed in previous cycle that will be the same if second half of previous polygon is the one the cell lies within
                if i2 - i:
                    normal2, j = compute_normal(i + step, relative_cell.GetX(), relative_cell.GetY())
                # auxiliary index i2, for normal2
                if step == -1:
                    i2 = i - 1
                else:
                    i2 = i + 2 * step
                # if neither of normals is less than allowed distance, new trajectory is inserted in between previous and following trajectory with help of normals ratio
                if np.round(normal2 / TSW):
                    # new alpha calculated from the normal ratio and new trajectory is generated (depends which direction the cycle is incrementing)
                    if step == -1:
                        new_alpha = TS[i][0] - abs(TS[i][0] - TS[i + step][0]) * normal2 / (normal1 + normal2)
                        TS.insert(i, [new_alpha, generate_trajectory(new_alpha)])
                    else:
                        new_alpha = TS[i][0] + abs(TS[i][0] - TS[i + step][0]) * normal1 / (normal1 + normal2)
                        TS.insert(i + step, [new_alpha, generate_trajectory(new_alpha)])
                    # index i needs to be set one less to start again at the same trajectory (or 2 more if incrementing from the end of set)
                    if step == -1:
                        i -= 2 * step
                    else:
                        i -= step
                    continue
                # if the following trajectory is the intersecting one, trajectory index is set accordingly
                else:
                    i += step
                    break
            break
    return trajectory_terrain_comparison(i, j, relative_cell)

def compute_normal(i, X_relative_cell, Y_relative_cell):
    """Computes perpendicular distance from closest segment of given trajectory and returns its size as well as index
    of first point of closest segment."""
    # Trajectory Point - Cell Distance List
    TPCDL = [((TS[i][1][0][j] - X_relative_cell) ** 2 + (TS[i][1][1][j] - Y_relative_cell) ** 2) ** (1 / 2) for j in range(len(TS[i][1][0]))]
    # index of closest point to cell is found
    j = TPCDL.index(min(TPCDL))
    # index of first of two closest points to cell is found
    if j != 0 and j != len(TPCDL) - 1:
        if TPCDL[j - 1] < TPCDL[j + 1]:
            j -= 1
    elif j == len(TPCDL) - 1:
        j -= 1
    # if closest point to cell mid point is closer than allowed distance, this distance is returned as there is no need to compute perpendicular distance to whole segment which can be only smaller than the point-cell distance
    if not np.round(min(TPCDL) / TSW):
        return min(TPCDL), j
    # calculate perpendicular distance (normal) from closest trajectory segment
    a = TPCDL[j]
    b = TPCDL[j + 1]
    c = ((TS[i][1][0][j] - TS[i][1][0][j + 1]) ** 2 + (TS[i][1][1][j] - TS[i][1][1][j + 1]) ** 2) ** (1 / 2)
    s = (a + b + c) / 2
    area = (s * (s - a) * (s - b) * (s - c)) ** (1 / 2)
    return area / c * 2, j

def trajectory_terrain_comparison(i, j, relative_cell):
    """Computes coordinates of terrain corresponding to each trajectory point and returns True or False depending
    on the result of terrain and trajectory point heights comparison."""
    # calculate azimuth of trajectory (shooting point to cell point)
    Azimuth = np.arctan((relative_cell.GetX() - SP.GetX())/(relative_cell.GetY() - SP.GetY()))
    # cycle iterates from first point of trajectory to the first point of segment closest to the cell point
    for X, Y in zip(TS[i][1][0][:j+1],TS[i][1][1][:j+1]):
        X_compare_point = SP.GetX() + X * np.sin(Azimuth)
        Y_compare_point = SP.GetY() + X * np.cos(Azimuth)
        Z_compare_point = int_function(X_compare_point, Y_compare_point)
        # if trajectory point is on or below terrain, False is returned
        if Y <= Z_compare_point:
            # there is a slight chance of the last point being within allowed distance radius of cell point which would make the cell hittable, not the opposite
            if not np.round(((X-relative_cell.GetX())**2+(Y-relative_cell.GetY())**2)**(1/2)/TSW):
                return True
            return False
    # if no trajectory point gets on or below terrain, cell is considered reachable, therefore True is returned
    return True

def plot_trajectory():
    import matplotlib.pyplot as plt  # na vykreslenie grafov
    plt.figure(figsize=(32, 18))
    for i in range(len(TS)):
        # plotting the points
        #plt.plot(TS[i][1][0], TS[i][1][1], markersize=5, linewidth=1, label=TS[i][0]/np.pi*180)
        plt.plot(TS[i][1][0], TS[i][1][1], '-', linewidth=1)

    #plt.plot(envelope[0], envelope[1], '-', linewidth=1)
    plt.plot([0, max(TS, key=lambda x: x[1][0][-1])[1][0][-1]], [DMAXH, DMAXH], '-', linewidth=1)

    end_cell = ogr.Geometry(ogr.wkbPoint)
    end_cell.AddPoint(SP.GetX(),SP.GetY()-100)
    profile = get_profile(end_cell)
    plt.plot(profile[0], profile[1], '-', linewidth=3)

    # plt.plot(ATF[0], ATF[1], '-', linewidth=3)
    # plt.plot(DTF[0], DTF[1], '-', linewidth=2)
    print(i)
    #plt.plot(temp_xyp[0], temp_xyp[1], 'r.', markersize=2)
    # xpar = max(TS, key=lambda x: x[1][0][-1])[1][0][-1]
    # ypar = max(TS[-1][1][1])
    # a = -ypar/xpar**2
    # b = 0
    # c = ypar
    # x = []
    # y = []
    # for i in range(0,int(xpar)+1):
    #     x.append(i)
    #     y.append(a*i**2+b*i+c)
    # plt.plot(x, y, 'b*', markersize=1)

    # ohranicenie, popis osi a nastavenie rovnakej mierky v smere oboch osi
    plt.xlim(0, profile[0][-1])
    plt.ylim(min(profile[1])-10, max(profile[1])+20)
    plt.xlabel("vzdialenosť [m]")
    plt.ylabel("výška [m]")
    plt.gca().set_aspect('equal', adjustable='box')

    # function to show the plot
    #plt.legend()
    plt.savefig('filename.png', dpi=900)
    plt.show()

def get_profile(end_cell):
    Azimuth = np.arctan((end_cell.GetX() - SP.GetX())/(end_cell.GetY() - SP.GetY()))
    profile = [[],[]]
    s = 0
    cell_dist = ((SP.GetX() - end_cell.GetX()) ** 2 + (SP.GetY() - end_cell.GetY()) ** 2) ** (1 / 2)
    while s < cell_dist:
        X_compare_point = SP.GetX() + s * np.sin(Azimuth)
        Y_compare_point = SP.GetY() + s * np.cos(Azimuth)
        Z_compare_point = int_function(X_compare_point, Y_compare_point)
        profile[0].append(s)
        profile[1].append(Z_compare_point)
        s += 0.5
    return profile



#######################################################################
## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr_clip2.tif" #path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\points1.shp"   #path to point layer
line_layer_path = None #r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp" #path to line layer, if obstacles are not to be used, set to None
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed"  #path to folder, where the file will be saved
throwshed_file = r"throwshedtestnew"   #name of output throwshed file

## SETTINGS
use_viewshed = 1 #utilization of viewshed, that will clip throwshed, No = 0, Yes = 1
band_number = 1 #selected band from DEM, default = 1
interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places)
EPSG = 8353 #EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG

## VARIABLES
initial_height = 1.7 #initial height of projectile above DEM when shot [m]
alpha_min = 0.0 #minimum of vertical angle range at which the projectile is shot [°]
alpha_max = 90.0 #maximum of vertical angle range at which the projectile is shot [°]
gravitational_acceleration = -9.81 #gravitational acceleration [m/s^2]
initial_velocity = 30 #initial velocity of projectile when shot [m/s]
air_density = 1.225 #air density [kg/m^3]
drag_coefficient = 2.0 #aerodynamic drag coefficient of projectile
cross_sectional_area = 0.000050 #cross-sectional area of the projectile [m^2]
mass = 0.035 #projectile mass [kg]
dalpha = 5 #step in vertical angle range [°]
trajectory_segment_width = None #distance step, at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step
eyes_height = 1.6 #shooter eye height above DEM for viewshed [m]
target_height = 1.7 #target height for viewshed [m]
wall_height = 4.0 #obstacle/wall height (if obstacle option is used) [m]
constant = 1 #constant multipling the drag coefficient within wobble distance of an arrow
area_addition = 0.0 #average addition to cross-sectional area of an arrow within wobble distance [m^2]
wobble_distance = 40 #wobble distance - distance at which an arrow stops wobbling [m]


main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, use_viewshed, EPSG,
         cumulative_throwshed, initial_height, initial_velocity, drag_coefficient, cross_sectional_area, mass,
         eyes_height, target_height, wall_height, constant, area_addition, wobble_distance, band_number=band_number,
         interpolation=interpolation, alpha_min=alpha_min, alpha_max=alpha_max,
         gravitational_acceleration=gravitational_acceleration, air_density=air_density, dalpha=dalpha,
         trajectory_segment_width=None)

# global IH, IV, DC, CSA, M, CONST, AA, WD, GA, AD, point_height, TSW, DMINH, DMAXH, AL
# IH, IV, DC, CSA, M, CONST, AA, WD, GA, AD, point_height, TSW, DMINH, DMAXH = initial_height, initial_velocity, drag_coefficient, \
#         cross_sectional_area, mass, constant, area_addition, wobble_distance, gravitational_acceleration, air_density, \
#                                                                0.0, 1, 0.0, 200.0

# throwshed()

#plot_trajectory()
#print(np.degrees(TS[iti-1][0]),np.degrees(TS[iti][0]),np.degrees(TS[iti+1][0]),np.degrees(TS[iti+2][0]),np.degrees(TS[-1][0]))

# for x1, y1, x2, y2 in zip(TS[iti][1][0],TS[iti][1][1],TS[iti+1][1][0],TS[iti+1][1][1]):
#     print(x1,y1,x2,y2)
#     if x2 > x1:
#         print('yes')


# import matplotlib.pyplot as plt  # na vykreslenie grafov
# plt.figure(figsize=(32, 18))
# x1 = np.multiply(TS[iti][1][0],1000000000)
# y1 = np.multiply(TS[iti][1][1],1000000000)
# x2 = np.multiply(TS[iti+1][1][0],1000000000)
# y2 = np.multiply(TS[iti+1][1][1],1000000000)
# # plotting the points
# plt.plot(x1, y1, '.-', markersize=4, linewidth=1, label='iti')
# plt.plot(x2, y2, '.-', markersize=4, linewidth=1, label='iti+1')
#
# # ohranicenie, popis osi a nastavenie rovnakej mierky v smere oboch osi
# xlim1 = (round(326.3099787361053,9)-0.000000005)*1000000000
# xlim2 = (round(326.3099787361053,9)+0.000000005)*1000000000
# ylim1 = (round(105.41666740186065,9)-0.000000005)*1000000000
# ylim2 = (round(105.41666740186065,9)+0.000000005)*1000000000
# plt.xlim(xlim1, xlim2)
# plt.ylim(ylim1, ylim2)
# plt.xlabel("vzdialenosť [m]")
# plt.ylabel("výška [m]")
# plt.gca().set_aspect('equal', adjustable='box')
#
# # function to show the plot
# plt.legend()
# plt.savefig('filename.png', dpi=900)
# plt.show()
#
# print(find_intersection(326.03612131613994,105.83548252333263,326.30997873207264,105.41666740802577,
#                         326.30997873610530,105.41666740186065,326.03612132016600,105.83548251716319))
#
# print(np.degrees(np.arctan((10541666740802577-10583548252333263)/(32630997873207264-32603612131613994))))
# print(np.degrees(np.arctan((10541666740186065-10583548251716319)/(32630997873610530-32603612132016600))))

#plot_trajectory()





"""
    # definicia vektorov, do ktorych sa budu ukladat suradnice najvzdialenejsich bodov jednotlivych azimutov
    X_coor_point_polygon = []
    Y_coor_point_polygon = []

    #Vytvorenie listu so vsetkymi hodnotami azimutov, ktore sa pouziju v cykle
    azimuth_list = np.arange(min_azimuth, max_azimuth + dazimuth, dazimuth)

    # Otacame pod azimutom (cyklus) a porovnavame hodnoty z y_r s DMR
    for Azimuth in azimuth_list:
        S = []  #vektor vzdialenosti k najvzdialenejsim bodom pri jednotlivych uhloch vystrelu
        # cyklus kde sa prestriedaju vsetky trajektorie (vsetky uhly vystrelu)
        for i in range(0,len(alfa_list)):
            j = 0
            r = 0
            #cyklus, kde sa meni vzdialenost
            while True:
                r += dr
                #vypocet suradnic so vzdialenostou dr-nasobku a pod Azimutom
                X_coor_compare_point = X_coor_point + r*np.sin(Azimuth)
                Y_coor_compare_point = Y_coor_point + r*np.cos(Azimuth)

                # Interpolacia DMR v porovnavanom bode
                if int_compare == 0:
                    #interpolacia DMR v bode porovnania (nearest neighbour)
                    dem_int_cell_column = round(np.abs((X_coor_compare_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1]))
                    dem_int_cell_row = round(np.abs((Y_coor_compare_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5]))
                    dem_int_point_height = dem_array[dem_int_cell_row][dem_int_cell_column]
                elif int_compare == 1:
                    #interpolacia DMR v bode porovnania (linear)
                    dem_int_cell_column = np.floor(np.abs((X_coor_compare_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1])).astype(np.int32) #najblizsii nizsi stlpec v array od bodu, X coor of the left column in set of four cells
                    dem_int_cell_row = np.floor(np.abs((Y_coor_compare_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5])).astype(np.int32)    #najblizsii nizsi riadok v array od bodu, Y coor of the upper row in set of four cells
                    X_coor_cell_1 = dem_gt[0] + dem_gt[1]/2 + dem_int_cell_column*dem_gt[1] #X suradnica stredov lavych buniek
                    # X_coor_cell_2 = dem_gt[0] + dem_gt[1]/2 + (dem_int_cell_column+1)*dem_gt[1] #X suradnica stredov pravych buniek
                    Y_coor_cell_1 = dem_gt[3] + dem_gt[5]/2 + (dem_int_cell_row+1)*dem_gt[5] #Y suradnica stredov dolnych buniek
                    # Y_coor_cell_2 = dem_gt[3] + dem_gt[5]/2 + dem_int_cell_row*dem_gt[5] #Y suradnica stredov hornych buniek
                    H_1 = dem_array[dem_int_cell_row][dem_int_cell_column]  #H lavej hornej bunky
                    H_2 = dem_array[dem_int_cell_row][dem_int_cell_column+1]  #H pravej hornej bunky
                    H_3 = dem_array[dem_int_cell_row+1][dem_int_cell_column]  #H lavej dolnej bunky
                    H_4 = dem_array[dem_int_cell_row+1][dem_int_cell_column+1]  #H pravej dolnej bunky
                    H_int_1 = ((X_coor_compare_point-X_coor_cell_1)*(H_4-H_3))/(np.abs(dem_gt[1])) + H_3   #Interpolovana vyska na dolnej linii
                    H_int_2 = ((X_coor_compare_point-X_coor_cell_1)*(H_2-H_1))/(np.abs(dem_gt[1])) + H_1   #Interpolovana vyska na hornej linii
                    dem_int_point_height = ((Y_coor_compare_point-Y_coor_cell_1)*(H_int_2-H_int_1))/(np.abs(dem_gt[5])) + H_int_1   #Interpolovana vyska medzi dolnou a hornou liniou
                else:
                    print("Hodnota int_compare neznama.")
                    exit()

                #porovnanie vysky bunky s vyskou sipu, ak je sip pod DMR, pripise sa maximalna vzdialenost pre konkretny uhol vystrelu
                if dem_int_point_height >= y_r[i][j]:
                    S.append(r)
                    break
                j += 1

            #nakoniec sa vyhlada maximalna vzdialenost spomedzi vsetkych pocitanych pre kazdy uhol vystrelu a zapisu sa suradnice najvzdialenejseho bodu pre dany azimut
            if i == range(0,len(alfa_list))[-1]:
                max_r = max(S)
                X_coor_point_polygon.append(X_coor_point + max_r*np.sin(Azimuth))
                Y_coor_point_polygon.append(Y_coor_point + max_r*np.cos(Azimuth))

    #######################################################################
    ## VYTVORENIE VYSTUPNEJ VRSTVY
    
    # vytvorenie novej geometrie
    throwshed_ring = ogr.Geometry(ogr.wkbLinearRing)
    # ak je azimut v celom rozsahu (throwshed pre cele okolie), pridaju sa len body dopadu
    if max_azimuth - min_azimuth == np.pi*2:
        for X, Y in zip(X_coor_point_polygon,Y_coor_point_polygon):
            throwshed_ring.AddPoint(X, Y)
    # ak je azimut iba v konkretnom rozsahu a nepocita sa throwshed pre cele okolie, treba pridat aj bod vystrelu, aby sa to spravne vykreslilo        
    else:
        throwshed_ring.AddPoint(X_coor_point, Y_coor_point)   #prvy bod (vystrelu) totozny s poslednym
        # pridanie zvysnych bodov do geometrie
        for X, Y in zip(X_coor_point_polygon,Y_coor_point_polygon):
            throwshed_ring.AddPoint(X, Y)
        throwshed_ring.AddPoint(X_coor_point, Y_coor_point)   #posledny bod (vystrelu) totozny s prvym

    # vytvorenie polygonu
    throwshed_polygon = ogr.Geometry(ogr.wkbPolygon)
    throwshed_polygon.AddGeometry(throwshed_ring)

    # ulozenie polygonu do vrstvy
    driver = ogr.GetDriverByName("ESRI Shapefile")
    throwshed_outds = driver.CreateDataSource(throwshed_output_folder + "\\" + throwshed_file + "_temp.shp")
    
    # definicia referencneho systemu
    srs = osr.SpatialReference()
    if keep_point_crs == 0:
        srs.ImportFromEPSG(EPSG)
    else:
        srs = point_layer.GetSpatialRef()
    throwshed_outlayer = throwshed_outds.CreateLayer(throwshed_file + "_temp", srs)
    
    # pridanie polygonu do feature a jeho ulozenie do vystupnej vrstvy
    throwshed_feature = ogr.Feature(throwshed_outlayer.GetLayerDefn())
    throwshed_feature.SetGeometry(throwshed_polygon)
    throwshed_outlayer.CreateFeature(throwshed_feature)

    # Vypocet maximalnej vzdialenosti pre viewshed z geoudajov vstupneho rastra a bodu vystrelu
    max_distance_4 =  (max([X_coor_point-dem_gt[0], dem_gt[0]+dem_gt[1]*len(dem_array[1])-X_coor_point])**2 + max([dem_gt[3]-Y_coor_point, Y_coor_point-(dem_gt[3]+dem_gt[5]*len(dem_array))])**2)**(1/2)
    
    # VYUZITIE VIEWSHED-U
    if use_viewshed == 1:
        # treba zavriet polygonovu vrstvu
        throwshed_outds = throwshed_outlayer = throwshed_feature = None
        # vytvorenie rastra viditelnosti, ulozi sa ako viewshed.tif do adresara s vystupnym throwshedom
        gdal.ViewshedGenerate(srcBand=dem_band, driverName='GTiff', targetRasterName=throwshed_output_folder + "\\" + viewshed_file + ".tif", creationOptions=None, observerX=X_coor_point, observerY=Y_coor_point, observerHeight=h_E+feet_height, targetHeight=h_T, visibleVal=1, invisibleVal=0, outOfRangeVal=0, noDataVal=-9999, dfCurvCoeff=0.85714, mode=2, maxDistance=np.max(max_distance_4))
        # otvorenie viewshed rastra
        viewshed_ds = gdal.Open(throwshed_output_folder + "\\" + viewshed_file + ".tif")
        # orezanie rastra viditelnosti throwshedom
        gdal.Warp(throwshed_output_folder + "\\" + throwshed_file + ".tif", viewshed_ds, cutlineDSName = throwshed_output_folder + "\\" + throwshed_file + "_temp" + ".shp", cropToCutline = False, dstNodata = 0)
        # vymazanie polygonu .shp s throwshedom a rastra .tif s viewshedom, tiez DMP sa vymaze
        for format in [file.split('.')[1] for file in os.listdir(throwshed_output_folder) if file.split('.')[0] == throwshed_file + "_temp"]:
            os.remove(throwshed_output_folder + '\\' + throwshed_file +"_temp." + format)
        viewshed_ds = None
        os.remove(throwshed_output_folder + "\\" + viewshed_file + ".tif")
    # NEVYUZITIE VIEWSHED-U
    elif use_viewshed == 0:
        # najprv treba vytvorit raster a dat mu nastavenia
        throwshed_ds = gdal.GetDriverByName('GTiff').Create(throwshed_output_folder + "\\" + throwshed_file + ".tif", xsize = dem_array.shape[1], ysize = dem_array.shape[0], bands = 1, eType = gdal.GDT_Byte)
        throwshed_ds.SetGeoTransform(dem_gt)
        throwshed_ds.SetProjection(srs.ExportToWkt())   #SS bude nastaveny ako bol aj pri polygonovej vrstve
        throwshed_band = throwshed_ds.GetRasterBand(1)
        throwshed_band.SetNoDataValue(0)
        # throwshed polygon sa rasterizuje, [1] - priradenie hodnot do pasma 1, burn_values=[1] - priradenie hodnot buniek = 1
        gdal.RasterizeLayer(throwshed_ds, [1], throwshed_outlayer, burn_values=[1])
        # nakoniec novovytvorena vrstva, datasource aj prvok treba dat rovne None, lebo inak sa nezobrazi spravne v QGISe
        throwshed_outds = throwshed_ds = throwshed_outlayer = throwshed_feature = None
        # vektorova podoba sa vymaze a zostane len rastrova, tiez DMP sa vymaze
        for format in [file.split('.')[1] for file in os.listdir(throwshed_output_folder) if file.split('.')[0] == throwshed_file + "_temp"]:
            os.remove(throwshed_output_folder + '\\' + throwshed_file +"_temp." + format)

    
    # SCITAVANIE RASTROV THROWSHED-OV VIACERYCH BODOV
    # prvy raster sa nacita do array ako zakladny
    if point_number_once == 0 and point_count > 1:
        # import prveho throwshedu
        throwshed_ds = gdal.Open(throwshed_output_folder + "\\" + throwshed_file + ".tif")
        # vyber pasma
        throwshed_band = throwshed_ds.GetRasterBand(1)
        # pridelenie hodnot rastra prveho throwshedu do hlavneho array
        throwshed_main_array = throwshed_band.ReadAsArray()
        point_number_once += 1  #zaruci ze sa zakladny raster do array ulozi iba raz
        throwshed_ds = throwshed_band = None
        # ak je bodov viac, a rastrov sa teda vytvori viac, je (zrejme) treba prvotny raster throwshedu vymazat
        os.remove(throwshed_output_folder + "\\" + throwshed_file + ".tif")

    # ak je bodov viac ako 1, vysledne rastre sa budu nacitavat a scitavat
    if point_number > 0:
        # import dalsieho throwshedu
        throwshed_ds = gdal.Open(throwshed_output_folder + "\\" + throwshed_file + ".tif")
        # vyber pasma
        throwshed_band = throwshed_ds.GetRasterBand(1)
        # pridelenie hodnot rastra dalsieho throwshedu do array
        throwshed_array = throwshed_band.ReadAsArray()
        throwshed_ds = throwshed_band = None
        # priebezne mazanie rastrov throwshedu jednotlivych bodov
        os.remove(throwshed_output_folder + "\\" + throwshed_file + ".tif")
        # pripocitanie array dalsieho rastra do hlavneho array
        throwshed_main_array = np.add(throwshed_main_array,throwshed_array)
        # ak bude cyklus riesit throwshed posledneho bodu, hlavny array sa ulozi do vysledneho rastra throwshedu
        if point_number == range(0,point_count)[-1]:
            # vytvorenie a nastavenie rastra vysledneho throwshedu
            throwshed_driver = gdal.GetDriverByName("GTiff")
            throwshed_outds = throwshed_driver.Create(throwshed_output_folder + "\\" + throwshed_file + ".tif", xsize = throwshed_main_array.shape[1], ysize = throwshed_main_array.shape[0], bands = 1, eType = gdal.GDT_Float32)
            throwshed_outds.SetGeoTransform(dem_gt)
            throwshed_outds.SetProjection(srs.ExportToWkt())
            throwshed_outband = throwshed_outds.GetRasterBand(1)
            # pridelenie kumulativneho alebo jednoducheho throwshedu
            if cumulative_throwshed == 0:
                throwshed_main_array = np.where((throwshed_main_array >= 1),1,0)
                throwshed_outband.WriteArray(throwshed_main_array)  #hodnoty rastra mozu byt 0 a 1
            elif cumulative_throwshed == 1:
                throwshed_outband.WriteArray(throwshed_main_array)  #hodnoty rastra mozu byt 0, 1, 2, 3...
            else:
                print("Volba kumulativneho throwshedu zle nastavena.")
                exit()
            throwshed_outband.SetNoDataValue(0)
            throwshed_outds = throwshed_outband = None

# zavretie a vymazanie DMP
if use_line == 1:
    dmp_outds = dem_band = None
    os.remove(throwshed_output_folder + "\\" + buffer_file + "_dmp_temp.tif")
"""