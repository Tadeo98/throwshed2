#######################################################################
## THROWSHED ##
#######################################################################
# Main #
#######################################################################

"""
Main throwshed analysis script that creates geotiff containing throwshed result.
"""

import os
import numpy as np
from osgeo import gdal, ogr, osr
from timeit import default_timer as timer
import numerical_methods

#######################################################################
## FUNCTIONS

def main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, throwshed_mode,
         use_viewshed, use_lines, cumulative_throwshed, EPSG, atmosphere_type, numerical_method,
         trajectory_segment_dimension, irregular_projectile, initial_height, initial_velocity, drag_to_mach, temperature, diameter, mass, cross_sectional_area,
         eyes_height, target_height, wall_height, wall_width, peak_drag, peak_area, oscillation_distance, oscillation_frequency, band_number=1,
         interpolation=1, alpha_min=-90.0, alpha_max=90.0, gravitational_acceleration=9.81, air_density=1.225,
         trajectory_segment_size=None):
    """Just main function with controls, global variables settings and triggers to other functions"""
    # making sure the vertical angle has correct range
    if alpha_max < alpha_min:
        print("Minimal vertical shooting angle higher than maximal.")
        exit()
    if not -90 <= alpha_max <= 90 or not -90 <= alpha_min <= 90:
        print("Minimal or maximal vertical shooting angle out of allowed range <-90°,+90°>.")
        exit()
    # Global variables
    global SRS, DP, PLP, TOF, TF, TM, UV, UL, CT, ATM, TSD, IP, NM, BN, INT, BF, TSS, RR, AL, DDS, DB, DA, DGT, DMINH, \
        DMAXH, IH, IV, D2M, T0, DIA, M, CSA, GA, AD, EH, TH, PD, PA, OD, OF, NDV, TA, VDS, VB, VA, VGT
    # CRS and other variable definition
    SRS = osr.SpatialReference()
    SRS.ImportFromEPSG(EPSG)
    TOF, TF, TM, UV, UL, CT, ATM, TSD, IP, BN, INT, IH, IV, D2M, T0, DIA, M, CSA, GA, AD, EH, TH, PD, PA, OD, OF = \
        throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed, use_lines, cumulative_throwshed, \
        atmosphere_type, trajectory_segment_dimension, irregular_projectile,\
        band_number, interpolation, initial_height, initial_velocity, drag_to_mach, temperature, diameter, mass, cross_sectional_area,\
        gravitational_acceleration, air_density, eyes_height, target_height, peak_drag, peak_area, oscillation_distance, oscillation_frequency
    # specific function that will compute ballistic trajectory is taken from numerical_methods module
    NM = getattr(numerical_methods, numerical_method)
    # get DEM data and assign them as global (and referencing datasource)
    DDS, DB, DA, DGT, NDV = get_raster_from_file(dem_path)
    # assign trajectory segment size
    TSS = np.min(np.abs([DGT[1],DGT[5]])) if trajectory_segment_size == None else trajectory_segment_size
    # raster resolution (cell size)
    RR = np.min(np.abs([DGT[1],DGT[5]]))
    # trajectory segment size cannot be larger than the size of a raster cell, and if so, will be set to raster cell size
    if TSS > RR:
        TSS = RR
    # obtain list of point geometries (and all referencing data it's dependent on)
    point_layer_ds, point_layer, point_feature_list, point_geom_list = get_geom_list_from_file(point_layer_path)
    # burn lines as obstacles into DEM, creating new DEM (Digital terrain model -> Digital surface model)
    if UL:
        burn_obstacles(line_layer_path, wall_height, wall_width)
    # get minimum and maximum DEM height
    DMINH, DMAXH = get_min_max_height()
    # obtain list of vertical angles
    AL = np.linspace(np.radians(alpha_min), np.radians(alpha_max), 37)
    # throwshed array containing zeroes at first, will be edited later, dimensions same as dimensions of DEM raster
    TA = [np.zeros((DA.shape[0], DA.shape[1]), np.int16), np.zeros((DA.shape[0], DA.shape[1]), np.int16)]
    # cycle calculating throwshed for each point
    for i, point_geom in enumerate(point_geom_list):
        # compute throwshed for 1 point
        throwshed(point_geom, i)
        # temporary viewshed files have to be removed
        if UV:
            VDS = VB = VA = VGT = None
            os.remove(TOF + "\\viewshed.tif")
    # finally, array is written into band of output throwshed raster
    create_raster_file(TF, TA, 1, gdal.GDT_Int16, NDV)
    DDS = DB = DA = DGT = NDV = None
    # Digital surface model file will be removed as well
    if UL:
        remove_temp_files(DSM)

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

def burn_obstacles(line_layer_path, wall_height, wall_width):
    """Lines that resemble walls/obstacles are burnt into DEM"""
    # get list of line geometries
    line_ds, line_layer, line_feature_list, line_geom_list = get_geom_list_from_file(line_layer_path)
    # calculate minimal buffer distance, at which the obstacle will always be respected (half of diagonal)
    min_buffer_dist = (DGT[1] ** 2 + DGT[5] ** 2) ** (1 / 2) / 2
    # if inserted wall width/2 is less than minimal buffer distance, buffer distance has to be set to this value
    buffer_dist = min_buffer_dist if wall_width/2 < min_buffer_dist else wall_width/2
    # 1st buffer is created and in case of multiple lines the cycle is utilized to create and unite buffers
    buffer_geom = line_geom_list[0].Buffer(buffer_dist, 1)
    if len(line_geom_list) > 1:
        for line_geom in line_geom_list[1:]:
            buffer_geom = buffer_geom.Union(line_geom.Buffer(buffer_dist, 1))
    # buffer file temporary name
    BFT = "buffer_temp"
    # create layer for buffer
    buffer_vector_ds, buffer_outlayer = create_and_use_outlayer(BFT, buffer_geom)
    # create buffer datasource
    buffer_raster_ds, buffer_band = create_raster_file(BFT, [0], 0, gdal.GDT_Float32, 0)
    # Buffer polygon is rasterized
    gdal.RasterizeLayer(buffer_raster_ds, [1], buffer_outlayer, burn_values=[wall_height])
    # Sum of initial dem and buffer rasters
    buffer_array = buffer_band.ReadAsArray()
    global DSM, DDS, DB, DA
    DA = np.add(DA, buffer_array)
    # create new raster datasource for DEM (DSM)
    DSM = "dsm_temp"
    DDS, DB = create_raster_file(DSM, [DA], 1, gdal.GDT_Float32, NDV)
    # delete all temporary files
    buffer_raster_ds = buffer_vector_ds = buffer_outlayer = buffer_band = None
    remove_temp_files(BFT)

def create_and_use_outlayer(layer_name, geom):
    """Creates file for output layer, which will contain new features. Returns created layer."""
    # create driver and output data source
    outds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(TOF + "\\" + layer_name + ".shp")
    # create output layer
    outlayer = outds.CreateLayer(layer_name, SRS)
    # feature definition and setting
    feature = ogr.Feature(outlayer.GetLayerDefn())
    feature.SetGeometry(geom)
    # assign feature into output layer
    outlayer.CreateFeature(feature)
    return outds, outlayer

def create_raster_file(raster_name, dem_array_list, method, GDT, no_data):
    """Creates raster file. Method 0 returns empty datasource and band. Method 1 returns datasource and band with written array"""
    # create driver and output data source
    outds = gdal.GetDriverByName('GTiff').Create(TOF + "\\" + raster_name + ".tif", xsize=DA.shape[1],
                                                 ysize=DA.shape[0], bands=len(dem_array_list), eType=GDT)
    # assign geotransformation, projection, band and nodata settings
    outds.SetGeoTransform(DGT)
    outds.SetProjection(SRS.ExportToWkt())
    for i, dem_array in enumerate(dem_array_list):
        raster_band = outds.GetRasterBand(i+1)
        raster_band.SetNoDataValue(no_data)
        if method:
            raster_band.WriteArray(dem_array)
    return outds, raster_band

def remove_temp_files(temp_file):
    """Deletes all temporary files with assigned name"""
    for format in [file.split('.',1)[1] for file in os.listdir(TOF) if file.split('.',1)[0] == temp_file]:
        os.remove(TOF + '\\' + temp_file + "." + format)

def throwshed(point_geom, k):
    """Calculates throwshed for 1 point"""
    global SP, WH
    # create shooting point with Z coordinate that is interpolated from DEM
    SP = ogr.Geometry(ogr.wkbPoint)
    SP.AddPoint(point_geom.GetX(), point_geom.GetY(), float(int_function(point_geom.GetX(), point_geom.GetY()))+IH)
    # generate set of trajectories for vertical angle range with basic step
    trajectory_simple_set()
    # insert trajectories between those from simple set, to ensure throwshed's edge accuracy
    trajectory_set()
    # define Ascending and Descending Trajectory Fields
    create_trajectory_fields()
    # if viewshed is ON, it has to be created and deleted later
    if UV:
        create_viewshed()
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
    # element in list begins with alpha value and continues with list of x and y (and z) coords lists in one trajectory, one of the numerical methods is chosen
    TS = [[alpha, NM(SP.GetZ(), AD, DIA, M, IV, alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)] for alpha in AL]

def trajectory_set():
    """Calculates and inserts trajectories between those in simple set, to make it denser and to ensure throwshed's
    edge accuracy. Calculates and returns trajectory envelope points list (its useful section)."""
    global TS, envelope
    # new and previous trajectory end x, first ones are random, just to make sure the cycle does not stop immediately
    ntex = [(max(TS, key=lambda x: x[1][0][-1])[1][0][-1]+RR)*2,(max(TS, key=lambda x: x[1][0][-1])[1][0][-1]+RR)*3]
    # cycle that finds the furthest possible trajectory for minimal DEM height respecting the edge accuracy
    while round(np.abs(ntex[0] - ntex[1])/RR):
        #most distant trajectory index
        mdti = TS.index((max(TS, key=lambda x: x[1][0][-1])))
        # adds new trajectories before and after current furthest trajectory
        if mdti != 0 and mdti != len(TS)-1:
            for new_alpha in [(TS[mdti+1][0] + TS[mdti][0]) / 2, (TS[mdti][0] + TS[mdti-1][0]) / 2]:
                TS.append([new_alpha, NM(SP.GetZ(), AD, DIA, M, IV, new_alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)])
        # for furthest trajectory that is also the first or last one, only one trajectory is added accordingly
        elif mdti == len(TS)-1:
            new_alpha = (TS[mdti][0] + TS[mdti-1][0]) / 2
            TS.append([new_alpha, NM(SP.GetZ(), AD, DIA, M, IV, new_alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)])
            mdti += 1 #this is just so that mdti gets higher like length of TS does (to cope with length of TS getting bigger, these 2 are compared in ntex)
        elif mdti == 0:
            new_alpha = (TS[mdti+1][0] + TS[mdti][0]) / 2
            TS.append([new_alpha, NM(SP.GetZ(), AD, DIA, M, IV, new_alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)])
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
    while True:
        # generate new trajectory in between actual and following
        new_alpha = (TS[iti][0] + TS[iti + 1][0]) / 2
        TS.insert(iti+1, [new_alpha, NM(SP.GetZ(), AD, DIA, M, IV, new_alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)])
        # in case of added trajectory having further reach than previously furthest trajectory (actual)
        if TS.index((max(TS, key=lambda x: x[1][0][-1]))) == iti + 1:
            iti += 1
            # even starting point of envelope needs to be updated
            envelope = [[TS[iti][1][0][-1]], [TS[iti][1][1][-1]]]
            XIIPR, YIIPR = envelope[0][0], envelope[1][0]
            continue

        # get intersection of actual and following (newly created) trajectory (Right Outer Intersection)
        XROI, YROI = calculate_intersection(TS[iti][1][0][1:], TS[iti][1][1][1:], TS[iti+1][1][0][1:], TS[iti+1][1][1][1:])
        # following intersections can be obtained only if the 2. following trajectory is not generated with 90-degrees shooting angle
        if TS[iti+2][0] != np.radians(90):
            # get intersection of following and 2. following trajectory (Left Outer Intersection)
            XLOI, YLOI = calculate_intersection(TS[iti+1][1][0][1:], TS[iti+1][1][1][1:], TS[iti+2][1][0][1:], TS[iti+2][1][1][1:])
            # get intersection of actual and 2. following trajectory (Inner Intersection), first one is calculated, the rest will be reused from outer intersections
            XII, YII = calculate_intersection(TS[iti][1][0][1:], TS[iti][1][1][1:], TS[iti+2][1][0][1:], TS[iti+2][1][1][1:])

        # if the last trajectory incorporated in cycle is the last one from the net with shooting angle value of 90 degrees.
        # This is because with 90 degrees trajectory left outer and inner intersections will be the same, which would lead to undesired behaviour
        if TS[iti+2][0] == np.radians(90):
            # if X coordinate of last point in newly created trajectory is less than cell size, last possible area of the net was made dense enough (used to be X of highest point, but last was chosen so that when searching for cell intersecting trajectories no trajectories are added if the cell falls between 2 last trajectories with the last having 90° angle, which would cause problems in rare situations)
            if not np.floor(TS[iti+1][1][0][-1]/RR):
                # update envelope with points from initial trajectory of last cycle, also update trajectory with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[iti].append(update_envelope(1, TS[iti], XIIPR, XROI, YROI, 0))
                # update last but one trajectory with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[iti+1].append(update_envelope(0, TS[iti+1], 0, 0, 0, YROI))
                # Highest Point Index of Last Trajectory
                HPILT = TS[-1][1][1].index((max(TS[-1][1][1])))
                envelope[0].append(0)
                envelope[1].append(TS[-1][1][1][HPILT])
                # update last trajectory (90 degrees one) with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[-1].append([[HPILT, HPILT], [len(envelope[0])-1, len(envelope[0])-1]])
                break
            # if not dense enough, density will be accomplished with new iteration
            else:
                continue

        # computation of Greatest Horizontal Width of Arc (triangle) which serves for stop criterion of densing iteration of one part of trajectory set
        if YROI < YII:
            # finds intersection of horizontal distance from inner intersection and particular segment on the arc of 1. following trajectory
            XI, YI = calculate_intersection([XII, XROI], [YII, YII], TS[iti+1][1][0], TS[iti+1][1][1])
            # Greatest Horizontal Width of Arc is calculated
            GHWA = XI - XII
        else:
            if YROI < YLOI:
                # finds intersection of horizontal distance from right intersection and particular segment on the arc of 2. following trajectory
                XI, YI = calculate_intersection([XLOI, XROI], [YROI, YROI], TS[iti+2][1][0], TS[iti+2][1][1])
            else:
                # finds intersection of horizontal distance from right intersection and particular segment on the arc of 1. following trajectory; average of XROI and XLOI because XROI itself would mean intersection at ROI
                XI, YI = calculate_intersection([XLOI, (XROI+XLOI)/2], [YROI, YROI], TS[iti+1][1][0], TS[iti+1][1][1])
            # Greatest Horizontal Width of Arc is calculated
            GHWA = XROI - XI

        # control whether arc (triangle) width (or height) criterion is met; when small enough, envelope is updated with new segments
        if round(GHWA / RR) and round((max(TS[iti+1][1][1]) - YII) / RR):
            continue
        else:
            # with each shooting point the amount of these inserted auxiliary trajectories would almost double which could create pointless amount of trajectories
            del TS[iti+1]
            # update envelope and update trajectory list with starting and ending point index of shared part between trajectory and envelope, indexes will be used when looking for cell neighbouring trajectories
            TS[iti].append(update_envelope(1, TS[iti], XIIPR, XII, YII, 0))
            # previous intersection for next cycle is assigned
            XIIPR, YIIPR = XII, YII
            # at least one of the conditions was met and the cycle can jump to next initial trajectory
            iti += 1
            # if the cycle comes to last trajectory, it breaks as there is no following trajectory
            if TS[iti][0] == AL[-1]:
                # even last trajectory needs to be updated with indexes and envelope will be finished
                TS[iti].append(update_envelope(0, TS[iti], 0, 0, 0, YII))
                break

def calculate_intersection(line1_x, line1_y, line2_x, line2_y):
    """Looks for intersection between two linestrings (trajectories or lines) based on coords of points in argument
    lists. Returns X and Y coordinates of the intersection."""
    # create 1. empty linestring geometry
    line1 = ogr.Geometry(ogr.wkbLineString)
    # fill the geometry with points
    for x, y in zip(line1_x, line1_y):
        line1.AddPoint(x, y)
    # create 2. empty linestring geometry
    line2 = ogr.Geometry(ogr.wkbLineString)
    # fill the 2. geometry with points
    for x, y in zip(line2_x, line2_y):
        line2.AddPoint(x, y)
    # compute intersection and return X and Y separately
    I = line1.Intersection(line2)
    return I.GetX(), I.GetY()

def update_envelope(method, TRAJ, XIIPR, XII, YII, YIPR):
    """Updates envelope with parts of trajectories and returns starting and ending indexes of points on part of
    trajectory that is shared with the envelope and also indexes of points on envelope that are boundaries of the
    shared trajectory part (previous and actual inner intersections). Trajectory points have to be between the envelope
    points."""
    global envelope
    # Reversed Trajectory list (X and Y coords)
    RT = [TRAJ[1][0][-1::-1], TRAJ[1][1][-1::-1]]
    # for regular parts
    if method:
        # on reversed trajectory, index of first point that will be incorporated within the envelope is found
        for i in range(len(RT[0]) - 1):
            if RT[0][i] >= XIIPR >= RT[0][i + 1]:
                # Envelope Part Starting Point Index
                EPSPI = i + 1
                break
        # on reversed trajectory, index of last point that will be incorporated within the envelope is found
        for i in range(EPSPI - 1, len(RT[0]) - 1):
            if RT[0][i] >= XII >= RT[0][i + 1]:
                # Envelope Part Ending Point Index
                EPEPI = i
                break
        # starting envelope index as the index of envelope's last point before update
        SEI = len(envelope[0]) - 1
        # condition for rare situation where both inner intersections could fall within one segment of trajectory
        if EPEPI >= EPSPI:
            # envelope is updated with all points between starting and ending point
            for i in range(EPSPI, EPEPI + 1):
                envelope[0].append(RT[0][i])
                envelope[1].append(RT[1][i])
        # ending envelope index as the index of envelope's last point after update
        EEI = len(envelope[0])
        # lastly, envelope is updated with the inner intersection point
        envelope[0].append(XII)
        envelope[1].append(YII)
        # return indexes of the first and last point of shared part, for trajectory direction of incrementing is from the left (shooting point), for envelope it's vice-versa
        return [[len(RT[0]) - 1 - EPSPI, len(RT[0]) - 1 - EPEPI], [SEI, EEI]]
    # for part of last but one trajectory when alpha of the last one is equal to 90° or part of last trajectory whose alpha is not equal to 90°
    else:
        # starting envelope index as the index of envelope's last point before update
        SEI = len(envelope[0]) - 1
        # starting and ending last but one trajectory point indexes
        EPSPI = 0
        # last points from last inserted trajectory and highest point of last trajectory are appended to envelope
        for i in range(len(RT[0])):
            if RT[1][i] > YIPR:
                if not EPSPI:
                    EPSPI = i
                envelope[0].append(RT[0][i])
                envelope[1].append(RT[1][i])
            if RT[1][i] == max(RT[1]):
                break
        EPEPI = i - 1
        # ending envelope index as the index of envelope's last point after update
        EEI = len(envelope[0]) - 1
        return [[len(RT[0]) - 1 - EPSPI, len(RT[0]) - 1 - EPEPI], [SEI, EEI]]

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
        # index of coords from last trajectory, where envelope connects with it
        i = TS[-1][2][0][1]
        ATF[0] = TS[0][1][0] + envelope[0] + TS[-1][1][0][i-1::-1]
        ATF[1] = TS[0][1][1] + envelope[1] + TS[-1][1][1][i-1::-1]
        # but for DTF with last trajectory shooting angle being 90 degrees, i is edited so that only last point down at min DEM height is added to polygon
        i = -2 if TS[-1][0] == np.radians(90) else i
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
        for j in range(DA.shape[1]):
            # with multiple shooting points nodata value can already be assigned to the cell, therefore the algorithm jumps to following cell
            if TA[0][i][j] == NDV:
                continue
            # nodata value is assigned to both arrays for both bands (ATF and DTF)
            if not k and DA[i][j] == NDV:
                TA[0][i][j] = TA[1][i][j] = NDV
                continue
            # for simple throwshed, if cell already has value 1, cycle continues with following cell, otherwise for cumulative throwshed, cell is assessed, descending trajectory does not need to be checked because if the cell is hittable by ascending trajectory, it has to be hittable by descending as well
            if k and not CT and TA[0][i][j]:
                continue
            # if viewshed is incorporated and particular cell is not visible, nothing is added to throwshed cell, and for visible cells the algorithm proceeds with assessment of cells
            if UV:
                # for case when only visible, hittable areas are sought
                if UV == 1 and not VA[i][j]:
                    continue
                # for case when only invisible, hittable areas are sought
                elif UV == -1 and VA[i][j]:
                    continue
            # calculate coordinates of cell's middle point and its horizontal distance from shooting point
            X_coor_cell = DGT[0]+(j+1/2)*DGT[1]
            Y_coor_cell = DGT[3]+(i+1/2)*DGT[5]
            cell_distance = ((SP.GetX() - X_coor_cell)**2 + (SP.GetY()-Y_coor_cell)**2)**(1/2)
            # create cell point with relative coordinates in the plane of trajectories and find out whether it's within the field, if so, further actions are conducted
            relative_cell = ogr.Geometry(ogr.wkbPoint)
            # also create cell point with absolute coordinates in the plane of projection plane, will be used in terrain comparison to calculate azimuth
            absolute_cell = ogr.Geometry(ogr.wkbPoint)
            relative_cell.AddPoint(cell_distance, float(DA[i][j]))
            absolute_cell.AddPoint(X_coor_cell, Y_coor_cell)
            # detect cell within the fields and call function to find cell intersecting trajectory and to determine whether the cell is reachable without any obstacles
            if ATF_polygon.Intersects(relative_cell):
                if TM:
                    if cell_availability(1, relative_cell, absolute_cell):
                        TA[0][i][j] += 1
                # for the case only cell's presence within the field is assessed
                else:
                    TA[0][i][j] += 1
            # can be None and also if CT is 0 and cell has a value already, no intersecting trajectory is found
            if DTF_polygon and not (not CT and TA[1][i][j]):
                if DTF_polygon.Intersects(relative_cell):
                    if TM:
                        if cell_availability(-1, relative_cell, absolute_cell):
                            TA[1][i][j] += 1
                    # for the case only cell's presence within the field is assessed
                    else:
                        TA[1][i][j] += 1

def cell_availability(dir, relative_cell, absolute_cell):
    """Finds trajectory that intersects the cell (or is close enough, within allowed distance). For ATF cycle
    increments from start to end of trajectory set and viceversa for DTF. Returns True if the cell is accessible
    or False if the cell is not accessible without any obstacles - this is determined by further function."""
    # Most Distant Trajectory Index
    MDTI = TS.index((max(TS, key=lambda x: x[1][0][-1])))
    # Zooming Index list containing indexes of assessed trajectories, for descending direction, first index is the one of the trajectory with furthest reach
    if dir == -1:
        ZI = [MDTI, int(MDTI + (len(TS)-MDTI) / 2), len(TS) - 1]
    else:
        ZI = [0, int(len(TS) / 2), len(TS) - 1]
    # cycle for zooming into the polygon of cell neighbouring trajectories
    while True:
        if dir == -1:
            # polygon also consists of envelope (starting and ending indexes of points of shared parts by trajectory and envelope are used)
            polygon = create_polygon_from_coords_list([envelope[0][TS[ZI[0]][2][1][0] + 1:TS[ZI[1]][2][1][0] + 1] +TS[ZI[1]][1][0][TS[ZI[1]][2][0][0] + 1:] +TS[ZI[0]][1][0][-1:TS[ZI[0]][2][0][0] - 1:-1], envelope[1][TS[ZI[0]][2][1][0] + 1:TS[ZI[1]][2][1][0] + 1] +TS[ZI[1]][1][1][TS[ZI[1]][2][0][0] + 1:] +TS[ZI[0]][1][1][-1:TS[ZI[0]][2][0][0] - 1:-1]])
        else:
            # very basic situation, envelope needs not to be used
            if ZI[1] <= MDTI:
                polygon = create_polygon_from_coords_list([TS[ZI[0]][1][0] + TS[ZI[1]][1][0][-1::-1], TS[ZI[0]][1][1] + TS[ZI[1]][1][1][-1::-1]])
            # situation where at least second trajectory is already intersecting other trajectories with further reach
            else:
                # situation where first of the trajectories is the one with furthest reach or the ones following
                if ZI[0] >= MDTI:
                    polygon = create_polygon_from_coords_list([TS[ZI[0]][1][0][:TS[ZI[0]][2][0][1]] + envelope[0][TS[ZI[0]][2][1][1]:TS[ZI[1]][2][1][1] + 1] +TS[ZI[1]][1][0][TS[ZI[1]][2][0][1] - 1::-1],TS[ZI[0]][1][1][:TS[ZI[0]][2][0][1]] + envelope[1][TS[ZI[0]][2][1][1]:TS[ZI[1]][2][1][1] + 1] +TS[ZI[1]][1][1][TS[ZI[1]][2][0][1] - 1::-1]])
                # situation where first of the trajectories precedes trajectory with furthest reach
                else:
                    polygon = create_polygon_from_coords_list([TS[ZI[0]][1][0] + envelope[0][:TS[ZI[1]][2][1][0] + 1] +TS[ZI[1]][1][0][TS[ZI[1]][2][0][0]::-1],TS[ZI[0]][1][1] + envelope[1][:TS[ZI[1]][2][1][0] + 1] +TS[ZI[1]][1][1][TS[ZI[1]][2][0][0]::-1]])
        if polygon.Intersects(relative_cell):
            ZI = [ZI[0], int((ZI[0] + ZI[1] ) / 2), ZI[1]]
        else:
            ZI = [ZI[1], int((ZI[1] + ZI[2]) / 2), ZI[2]]
        if abs(ZI[0] - ZI[2]) == 1:
            i = ZI[0]
            break

    # auxiliary indexes
    i1, i2 = -1, -1
    # Intersecting Trajectory Found list informing whether the previous and following intersecting trajectories were already found and compared with the terrain
    ITF = [1,1]
    # Inserted Trajectories Starting Index from which newly added trajectories will be removed right before returning from this function
    ITSI = -1
    # Inserted Trajectories Index Span across which the newly added trajectories will be removed
    ITIS = 0
    # at first, 2 surrounding trajectories are found by making polygon out of them and asking whether the cell lies within
    while True:
        # create polygon with existing/inserted trajectories
        polygon = create_polygon_from_coords_list([TS[i][1][0] + TS[i + 1][1][0][-1::-1], TS[i][1][1] + TS[i + 1][1][1][-1::-1]])
        # if the cell lies within, normals are computed to assess the smallest perpendicular distance from trajectory to cell
        if polygon.Intersects(relative_cell):
            # Inserted Trajectories Starting Index
            if ITSI == -1:
                ITSI = i+1
            # first condition is just to recycle normal1/2 computed in previous cycle that will be the same if first/second half of previous polygon is the one the cell lies within
            if i1 - i and ITF[0]:
                # compute normal from previous/following trajectory segment to cell
                normal1 = compute_normal(i, relative_cell.GetX(), relative_cell.GetY())
                # if cell is not close enough to trajectory to consider it as piercing trajectory, following trajectory is tested. If it is close enough, it will be used as the index of intersecting trajectory in the terrain comparison
                if not np.round(normal1 / RR):
                    if trajectory_terrain_comparison(i, relative_cell, absolute_cell):
                        del TS[ITSI:ITSI + ITIS]
                        return True
                    # 1 is changed to 0 so next time the condition is eluded
                    ITF[0] -= 1
            if i2 - i and ITF[1]:
                normal2 = compute_normal(i+1, relative_cell.GetX(), relative_cell.GetY())
                if not np.round(normal2 / RR):
                    if trajectory_terrain_comparison(i+1, relative_cell, absolute_cell):
                        del TS[ITSI:ITSI + ITIS]
                        return True
                    ITF[1] -= 1
            # if trajectories from both sides are intersecting and none of them returned True meaning reachable cell, cell is considered unreachable
            if not any(ITF):
                del TS[ITSI:ITSI+ITIS]
                return False
            # ratio for angle addition to angle of previous trajectory, if one of the trajectories is already intersecting, second one is being searched by halving the angle difference of surrounding trajectories, because by normal ratio new trajectory can fall on wrong side of the cell, to the one that has already been assessed, which will slow down the computation
            ratio = 1 / 2 if not ITF[0] or not ITF[1] else normal1 / (normal1 + normal2)
            # new alpha calculated from the ratio and new trajectory is generated
            new_alpha = TS[i][0] + (TS[i + 1][0] - TS[i][0]) * ratio
            TS.insert(i + 1, [new_alpha, NM(SP.GetZ(), AD, DIA, M, IV, new_alpha, T0, ATM, D2M, GA, TSS, TSD, IP, CSA, PD, PA, OD, OF, DMINH)])
            # index i needs to be set one less to start again at the same trajectory
            # auxiliary index i1/2 to find out if the normal1/2 was already computed, will be used in next iteration
            i -= 1
            i1 = i + 1
            i2 = i + 2
            # 1 trajectory added to the index span
            ITIS += 1
        i += 1

def compute_normal(i, X_relative_cell, Y_relative_cell):
    """Computes perpendicular distance from closest segment of given trajectory and returns its size as well as index
    of first point of closest segment."""
    # Trajectory Point - Cell Distance List
    TPCDL = [((TS[i][1][0][j] - X_relative_cell) ** 2 + (TS[i][1][1][j] - Y_relative_cell) ** 2) ** (1 / 2) for j in range(len(TS[i][1][0]))]
    # if closest point to cell mid point is closer than allowed distance, this distance is returned as there is no need to compute perpendicular distance to whole segment which can be only smaller than the point-cell distance
    if not np.round(min(TPCDL) / RR):
        return min(TPCDL)
    # index of closest point to cell is found
    j = TPCDL.index(min(TPCDL))
    # index of first of two closest points to cell is found
    if j != 0 and j != len(TPCDL) - 1:
        if TPCDL[j - 1] < TPCDL[j + 1]:
            j -= 1
    elif j == len(TPCDL) - 1:
        j -= 1
    # calculate perpendicular distance (normal) from closest trajectory segment
    a = TPCDL[j]
    b = TPCDL[j + 1]
    c = ((TS[i][1][0][j] - TS[i][1][0][j + 1]) ** 2 + (TS[i][1][1][j] - TS[i][1][1][j + 1]) ** 2) ** (1 / 2)
    s = (a + b + c) / 2
    area = (s * (s - a) * (s - b) * (s - c)) ** (1 / 2)
    return area / c * 2

def trajectory_terrain_comparison(i, relative_cell, absolute_cell):
    """Computes coordinates of terrain corresponding to each trajectory point and returns True or False depending
    on the result of terrain and trajectory point heights comparison."""
    # calculate azimuth of trajectory (shooting point to cell point), there is a chance of Y difference to be 0, therefore the exception
    dX = absolute_cell.GetX() - SP.GetX()
    dY = absolute_cell.GetY() - SP.GetY()
    # for case when the shooting point is right on teh cell, cell is automatically reachable
    if not round(dX/RR) and not round(dY/RR):
        return True
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
    # cycle iterates from first point of trajectory to the first point of segment closest to the cell point
    for X, Y in zip(TS[i][1][0],TS[i][1][1]):
        # if compared point is already above/below destination cell's area or beyond, cell is considered reachable
        if X >= relative_cell.GetX()-RR/2:
            return True
        X_compare_point = SP.GetX() + X * np.sin(Azimuth)
        Y_compare_point = SP.GetY() + X * np.cos(Azimuth)
        Z_compare_point = int_function(X_compare_point, Y_compare_point)
        # if trajectory point is on or below terrain, False is returned
        if Y <= Z_compare_point:
            return False

def create_viewshed():
    """Computes viewshed for one point which is saved temporarily, then load as an array."""
    global VDS, VB, VA, VGT
    # generate viewshed and save it as temporary file to throwshed directory
    gdal.ViewshedGenerate(srcBand=DB, driverName='GTiff', targetRasterName=TOF + "\\viewshed.tif", creationOptions=[], observerX=SP.GetX(), observerY=SP.GetY(), observerHeight=EH, targetHeight=TH, visibleVal=1, invisibleVal=0, outOfRangeVal=0, noDataVal=NDV, dfCurvCoeff=0.85714, mode=2, maxDistance=0)
    # open viewshed raster, Viewshed Array will be crucial
    VDS, VB, VA, VGT, ndv = get_raster_from_file(TOF + "\\viewshed.tif")
