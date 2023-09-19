#######################################################################
## THROWSHED ##
#######################################################################
import os
import numpy as np
from osgeo import gdal, ogr, osr
from timeit import default_timer as timer

#######################################################################
## FUNCTIONS

def main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed, use_lines, EPSG,
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
    global SRS, DP, PLP, TOF, TF, TM, UV, UL, CT, BN, INT, BF, TSW, AL, DDS, DB, DA, DGT, DMINH, DMAXH, IH, IV, \
        DC, CSA, M, CONST, AA, WD, GA, AD, EH, TH, NDV, TA, VDS, VB, VA, VGT
    # CRS and other variable definition
    SRS = osr.SpatialReference()
    SRS.ImportFromEPSG(EPSG)
    TOF, TF, TM, UV, UL, CT, BN, INT, IH, IV, DC, CSA, M, CONST, AA, WD, GA, AD, EH, TH = \
        throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed, use_lines, cumulative_throwshed, \
        band_number, interpolation, initial_height, initial_velocity, drag_coefficient, cross_sectional_area, mass, \
        constant, area_addition, wobble_distance, gravitational_acceleration, air_density, eyes_height, target_height
    # get DEM data and assign them as global (and referencing datasource)
    DDS, DB, DA, DGT, NDV = get_raster_from_file(dem_path)
    # assign trajectory segment width
    TSW = np.min(np.abs([DGT[1],DGT[5]])) if trajectory_segment_width == None else trajectory_segment_width
    # obtain list of point geometries (and all referencing data it's dependent on)
    point_layer_ds, point_layer, point_feature_list, point_geom_list = get_geom_list_from_file(point_layer_path)
    # burn lines as obstacles into DEM, creating new DEM (Digital terrain model -> Digital surface model)
    if UL:
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
        #break #throwshed only for the first point is computed
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

def burn_obstacles(line_layer_path, wall_height):
    """Lines that resemble walls/obstacles are burnt into DEM"""
    # get list of line geometries
    line_ds, line_layer, line_feature_list, line_geom_list = get_geom_list_from_file(line_layer_path)
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
    for format in [file.split('.')[1] for file in os.listdir(TOF) if file.split('.')[0] == temp_file]:
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
    trajectory_initial_set()
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
            # if X coordinate of last point in newly created trajectory is less than TSW, last possible area of the net was made dense enough (used to be X of highest point, but last was chosen so that when searching for cell intersecting trajectories no trajectories are added if the cell falls between 2 last trajectories with the last having 90° angle, which would cause problems in rare situations)
            if not np.floor(TS[iti+1][1][0][-1]/TSW):
                # Initial Trajectory Reversed list (X and Y coords)
                ITR = [TS[iti][1][0][-1::-1], TS[iti][1][1][-1::-1]]
                # update envelope with points from initial trajectory of last cycle, also update trajectory with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[iti].append(update_envelope(1, ITR, XIIPR, XROI, YROI, 0))
                # Last Inserted Trajectory Reversed list (X and Y coords)
                LITR = [TS[iti+1][1][0][-1::-1], TS[iti+1][1][1][-1::-1]]
                # update last but one trajectory with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[iti+1].append(update_envelope(0, LITR, 0, 0, 0, YROI))
                envelope[0].append(0)
                envelope[1].append(max(TS[-1][1][1]))
                # update last trajectory (90 degrees one) with shared part starting and ending point indexes (on trajectory as on envelope)
                TS[iti + 2].append([[TS[iti+2][1][1].index((max(TS[iti+2][1][1]))), TS[iti+2][1][1].index((max(TS[iti+2][1][1])))], [len(envelope[0])-1, len(envelope[0])-1]])
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
            # update envelope and update trajectory list with starting and ending point index of shared part between trajectory and envelope, indexes will be used when looking for cell neighbouring trajectories
            TS[iti].append(update_envelope(1, ITR, XIIPR, XII, YII, 0))
            # previous intersection for next cycle is assigned
            XIIPR, YIIPR = XII, YII
            # at least one of the conditions was met and the cycle can jump to next initial trajectory
            iti += 1
            # if the cycle comes to last trajectory, it breaks as there is no following trajectory
            if TS[iti][0] == AL[-1]:
                # even last trajectory needs to be updated with indexes
                ITR = [TS[iti][1][0][-1::-1], TS[iti][1][1][-1::-1]]
                TS[iti].append(update_envelope(0, ITR, 0, 0, 0, YII))
                break

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

def update_envelope(method, ITR, XIIPR, XII, YII, YROI):
    """Updates envelope with parts of trajectories and returns starting and ending indexes of points on part of
    trajectory that is shared with the envelope and also indexes of points on envelope that are boundaries of the
    shared trajectory part (previous and actual inner intersections). Trajectory points have to be between the envelope
    points."""
    global envelope
    # for regular parts
    if method:
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
        # starting envelope index as the index of envelope's last point before update
        SEI = len(envelope[0]) - 1
        # condition for rare situation where both inner intersections could fall within one segment of trajectory
        if EPEPI >= EPSPI:
            # envelope is updated with all points between starting and ending point
            for i in range(EPSPI, EPEPI + 1):
                envelope[0].append(ITR[0][i])
                envelope[1].append(ITR[1][i])
        # ending envelope index as the index of envelope's last point after update
        EEI = len(envelope[0])
        # lastly, envelope is updated with the inner intersection point
        envelope[0].append(XII)
        envelope[1].append(YII)
        # return indexes of the first and last point of shared part, for trajectory direction of incrementing is from the left (shooting point), for envelope it's vice-versa
        return [[len(ITR[0]) - 1 - EPSPI, len(ITR[0]) - 1 - EPEPI], [SEI, EEI]]
    # for part of last but one trajectory when alpha of the last one is equal to 90° or part of last trajectory whose alpha is not equal to 90°
    else:
        # starting envelope index as the index of envelope's last point before update, for last but one trajectory
        SEI = len(envelope[0]) - 1
        # starting and ending last but one trajectory point indexes
        EPSPI = 0
        # last points from last inserted trajectory and highest point of last trajectory are appended to envelope
        for i in range(len(ITR[0])):
            if ITR[1][i] > YROI:
                if not EPSPI:
                    EPSPI = i
                envelope[0].append(ITR[0][i])
                envelope[1].append(ITR[1][i])
            if ITR[0][i] == ITR[0][ITR[1].index((max(ITR[1])))]:
                break
        EPEPI = i - 1
        # ending envelope index as the index of envelope's last point after update, for last but one trajectory
        EEI = len(envelope[0]) - 1
        return [[len(ITR[0]) - 1 - EPSPI, len(ITR[0]) - 1 - EPEPI], [SEI, EEI]]

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



            start1 = timer()



            # with multiple shooting points nodata value can already be assigned to the cell, therefore the algorithm jumps to following cell
            if TA[0][i][j] == NDV:
                continue
            # nodata value is assigned to both arrays for both bands (ATF and DTF)
            if not k and DA[i][j] == NDV:
                TA[0][i][j] = TA[1][i][j] = NDV
                continue
            # for simple throwshed, if cell already has value 1, cycle continues with following cell, otherwise for cumulative throwshed, cell is assessed
            if k and CT == 0 and TA[0][i][j]:
                continue
            # if viewshed is incorporated and particular cell is not visible, nothing is added to throwshed cell, and for visible cells the algorithm proceeds with assessment of cells
            if UV:
                if not VA[i][j]:
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
                    if find_intersecting_trajectory(1, relative_cell, absolute_cell):
                        TA[0][i][j] += 1
                # for the case only cell's presence within the field is assessed
                else:
                    TA[0][i][j] += 1
            # can be None
            if DTF_polygon:
                if DTF_polygon.Intersects(relative_cell):
                    if TM:
                        if find_intersecting_trajectory(-1, relative_cell, absolute_cell):
                            TA[1][i][j] += 1
                    # for the case only cell's presence within the field is assessed
                    else:
                        TA[1][i][j] += 1

            end1 = timer()
            print(f'Cell {i} {j} took {end1-start1} seconds.')
            print()

def find_intersecting_trajectory(dir, relative_cell, absolute_cell):
    """Finds trajectory that intersects the cell (or is close enough, within allowed distance). For ATF cycle
    increments from start to end of trajectory set and viceversa for DTF. Returns True if the cell is accessible
    or False if the cell is not accessible without any obstacles - this is determined further function."""


    start4 = timer()


    # Most Distant Trajectory Index
    MDTI = TS.index((max(TS, key=lambda x: x[1][0][-1])))
    if dir == -1:
        # Zooming Index List containing indexes of assessed trajectories, for descending direction, first index is the one of the trajectory with furthest reach
        ZIL = [MDTI, int(MDTI + (len(TS)-MDTI) / 2), len(TS) - 1]
    else:
        ZIL = [0, int(len(TS) / 2), len(TS) - 1]
    # cycle for zooming into the polygon of cell neighbouring trajectories
    while True:
        for j in [0, 1]:
            if dir == -1:
                # polygon also consists of envelope (starting and ending indexes of points of shared parts by trajectory and envelope are used)
                polygon = create_polygon_from_coords_list([envelope[0][TS[ZIL[j]][2][1][0]+1:TS[ZIL[j+1]][2][1][0]+1] + TS[ZIL[j + 1]][1][0][TS[ZIL[j+1]][2][0][0]+1:] + TS[ZIL[j]][1][0][-1:TS[ZIL[j]][2][0][0]-1:-1], envelope[1][TS[ZIL[j]][2][1][0]+1:TS[ZIL[j+1]][2][1][0]+1] + TS[ZIL[j + 1]][1][1][TS[ZIL[j+1]][2][0][0]+1:] + TS[ZIL[j]][1][1][-1:TS[ZIL[j]][2][0][0]-1:-1]])
            else:
                # very basic situation, envelope needs not to be used
                if ZIL[j+1] <= MDTI:
                    polygon = create_polygon_from_coords_list([TS[ZIL[j]][1][0] + TS[ZIL[j + 1]][1][0][-1::-1], TS[ZIL[j]][1][1] + TS[ZIL[j + 1]][1][1][-1::-1]])
                # situation where at least second trajectory is already intersecting other trajectories with further reach
                else:
                    # situation where first of the trajectories is the one with furthest reach or the ones following
                    if ZIL[j] >= MDTI:
                        polygon = create_polygon_from_coords_list([TS[ZIL[j]][1][0][:TS[ZIL[j]][2][0][1]] + envelope[0][TS[ZIL[j]][2][1][1]:TS[ZIL[j+1]][2][1][1]+1] + TS[ZIL[j + 1]][1][0][TS[ZIL[j+1]][2][0][1]-1::-1], TS[ZIL[j]][1][1][:TS[ZIL[j]][2][0][1]] + envelope[1][TS[ZIL[j]][2][1][1]:TS[ZIL[j+1]][2][1][1]+1] + TS[ZIL[j + 1]][1][1][TS[ZIL[j+1]][2][0][1]-1::-1]])
                    # situation where first of the trajectories precedes trajectory with furthest reach
                    else:
                        polygon = create_polygon_from_coords_list([TS[ZIL[j]][1][0] + envelope[0][:TS[ZIL[j+1]][2][1][0]+1] + TS[ZIL[j + 1]][1][0][TS[ZIL[j+1]][2][0][0]::-1], TS[ZIL[j]][1][1] + envelope[1][:TS[ZIL[j+1]][2][1][0]+1] + TS[ZIL[j + 1]][1][1][TS[ZIL[j+1]][2][0][0]::-1]])
            if polygon.Intersects(relative_cell):
                break
        if abs(ZIL[j + 1] - ZIL[j]) == 1:
            i = ZIL[j]
            break
        ZIL = [ZIL[j], int(ZIL[j] + (ZIL[j + 1] - ZIL[j]) / 2), ZIL[j + 1]]


    end4 = timer()
    print(f'Zooming took {end4 - start4} seconds.')


    start5 = timer()

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
        # if the cell lies within, normals are computed to assess the smallest perpendicular distance from trajectory to cell
        if polygon.Intersects(relative_cell):
            # Inserted Trajectories Starting Index
            if ITSI == -1:
                ITSI = i+1
            # first condition is just to recycle normal1/2 computed in previous cycle that will be the same if first/second half of previous polygon is the one the cell lies within
            if i1 - i and ITF[0]:
                # compute normal from previous/following trajectory segment to cell
                normal1 = compute_normal(i, relative_cell.GetX(), relative_cell.GetY())
                # if cell is not close enough to trajectory to consider it as piercing trajectory, following trajectory is tested. If it is close enough, i will be used as the index of intersecting trajectory in the terrain comparison
                if not np.round(normal1 / TSW):
                    if trajectory_terrain_comparison(i, relative_cell, absolute_cell):
                        del TS[ITSI:ITSI + ITIS]


                        end5 = timer()
                        print(f'Intersecting took {end5 - start5} seconds.')


                        return True
                    # 1 is changed to 0 so next time the condition is eluded
                    ITF[0] -= 1
            if i2 - i and ITF[1]:
                normal2 = compute_normal(i+1, relative_cell.GetX(), relative_cell.GetY())
                if not np.round(normal2 / TSW):
                    if trajectory_terrain_comparison(i+1, relative_cell, absolute_cell):
                        del TS[ITSI:ITSI + ITIS]


                        end5 = timer()
                        print(f'Intersecting took {end5 - start5} seconds.')


                        return True
                    ITF[1] -= 1
            # if trajectories from both sides are intersecting and none of them returned True meaning reachable cell, cell is considered unreachable
            if not any(ITF):
                del TS[ITSI:ITSI+ITIS]


                end5 = timer()
                print(f'Intersecting took {end5 - start5} seconds.')


                return False
            # ratio for angle addition to angle of previous trajectory, if one of the trajectories is already intersecting, second one is being searched by halving the angle difference of surrounding trajectories, because by normal ratio new trajectory can fall on wrong side of the cell, to the one that has already been assessed, which will slow down the computation
            ratio = 1 / 2 if not ITF[0] or not ITF[1] else normal1 / (normal1 + normal2)
            # new alpha calculated from the ratio and new trajectory is generated
            new_alpha = TS[i][0] + (TS[i + 1][0] - TS[i][0]) * ratio
            TS.insert(i + 1, [new_alpha, generate_trajectory(new_alpha)])
            # index i needs to be set one less to start again at the same trajectory
            # auxiliary index i1/2 to find out if the normal1/2 was already computed, will be used in next iteration
            i -= 1
            i1 = i
            i2 = i + 2
            # 1 trajectory added to the index span
            ITIS += 1
        i += 1
        # create polygon with inserted trajectory/ies
        polygon = create_polygon_from_coords_list([TS[i][1][0] + TS[i + 1][1][0][-1::-1], TS[i][1][1] + TS[i + 1][1][1][-1::-1]])

def compute_normal(i, X_relative_cell, Y_relative_cell):
    """Computes perpendicular distance from closest segment of given trajectory and returns its size as well as index
    of first point of closest segment."""

    start2 = timer()


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




        end2 = timer()
        print(f'Normal took {end2 - start2} seconds.')


        return min(TPCDL)
    # calculate perpendicular distance (normal) from closest trajectory segment
    a = TPCDL[j]
    b = TPCDL[j + 1]
    c = ((TS[i][1][0][j] - TS[i][1][0][j + 1]) ** 2 + (TS[i][1][1][j] - TS[i][1][1][j + 1]) ** 2) ** (1 / 2)
    s = (a + b + c) / 2
    area = (s * (s - a) * (s - b) * (s - c)) ** (1 / 2)


    end2 = timer()
    print(f'Normal took {end2-start2} seconds.')

    return area / c * 2

def trajectory_terrain_comparison(i, relative_cell, absolute_cell):
    """Computes coordinates of terrain corresponding to each trajectory point and returns True or False depending
    on the result of terrain and trajectory point heights comparison."""


    start3 = timer()


    # calculate azimuth of trajectory (shooting point to cell point), there is a chance of Y difference to be 0, therefore the exception
    dX = absolute_cell.GetX() - SP.GetX()
    dY = absolute_cell.GetY() - SP.GetY()
    # for the case where shooting point and middle point of assessed cell are same, automatically reachable
    if not dX and not dY:
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
        if X >= relative_cell.GetX()-TSW/2:


            end3 = timer()
            print(f'Terrain comparison took {end3 - start3} seconds.')


            return True
        X_compare_point = SP.GetX() + X * np.sin(Azimuth)
        Y_compare_point = SP.GetY() + X * np.cos(Azimuth)
        Z_compare_point = int_function(X_compare_point, Y_compare_point)
        # if trajectory point is on or below terrain, False is returned
        if Y <= Z_compare_point:


            end3 = timer()
            print(f'Terrain comparison took {end3 - start3} seconds.')


            return False

def create_viewshed():
    """Computes viewshed for one point which is saved temporarily, then load as an array."""
    global VDS, VB, VA, VGT
    # generate viewshed and save it as temporary file to throwshed directory
    gdal.ViewshedGenerate(srcBand=DB, driverName='GTiff', targetRasterName=TOF + "\\viewshed.tif", creationOptions=[], observerX=SP.GetX(), observerY=SP.GetY(), observerHeight=EH, targetHeight=TH, visibleVal=1, invisibleVal=0, outOfRangeVal=0, noDataVal=NDV, dfCurvCoeff=0.85714, mode=2, maxDistance=0)
    # open viewshed raster, Viewshed Array will be crucial
    VDS, VB, VA, VGT, ndv = get_raster_from_file(TOF + "\\viewshed.tif")

def plot_trajectory(relative_cell, absolute_cell,jj,row,col,dir, poly_list, j, zoom):
    import matplotlib.pyplot as plt  # na vykreslenie grafov
    plt.figure(figsize=(32, 18))
    for i in range(len(TS)):
        # plotting the points
        #plt.plot(TS[i][1][0], TS[i][1][1], markersize=5, linewidth=1, label=TS[i][0]/np.pi*180)
        plt.plot(TS[i][1][0], TS[i][1][1], '.', markersize=1)

    #plt.plot(envelope[0], envelope[1], '-', linewidth=1)
    #plt.plot([0, max(TS, key=lambda x: x[1][0][-1])[1][0][-1]], [DMAXH, DMAXH], '-', linewidth=1)

    # end_cell = ogr.Geometry(ogr.wkbPoint)
    # end_cell.AddPoint(-488475.5,-1259010.5)
    profile = get_profile(absolute_cell)
    plt.plot(profile[0], profile[1], '-', linewidth=3)

    # plt.plot(TS[jj][1][0], TS[jj][1][1], '-', linewidth=1)
    # plt.plot(TS[jj+1][1][0], TS[jj+1][1][1], '-', linewidth=1)

    plt.plot(poly_list[0], poly_list[1], '-', linewidth=2)

    plt.plot(relative_cell.GetX(), relative_cell.GetY(), 'o', markersize=2)

    # plt.plot(ATF[0], ATF[1], '-', linewidth=2)
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
    plt.xlim(0, 155)
    plt.ylim(DMINH, max(TS[-1][1][1]))

    # plt.xlim(min(poly_list[0]), max(poly_list[0]))
    # plt.ylim(min(poly_list[1]), max(poly_list[1]))


    plt.xlabel("vzdialenosť [m]")
    plt.ylabel("výška [m]")
    plt.gca().set_aspect('equal', adjustable='box')

    # function to show the plot
    #plt.legend()


    plt.savefig(f'filename{row}_{col}_{dir}_zoom{zoom}_j{j}.png', dpi=300)

    start = timer()
    plt.show()
    end = timer()
    print('cas vykreslenia:', end-start)

def get_profile(end_cell):
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

#######################################################################
## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr_clip.tif" #path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\points1.shp"   #path to point layer
line_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp" #path to line layer
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed"  #path to folder, where the file will be saved
throwshed_file = r"throwshedtestnew"   #name of output throwshed file

## SETTINGS
throwshed_mode = 1 #what type of throwshed will be calculated, simple safety zone (cells within safety field) = 0, regular throwshed with trajectory assessment = 1
use_viewshed = 0 #utilization of viewshed, that will clip throwshed, No = 0, Yes = 1
use_lines = 0 #utilization of line layer, where lines serve as obstacles or walls and will be burnt into DEM, No = 0, Yes = 1
band_number = 1 #selected band from DEM, default = 1
interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places)
EPSG = 8353 #EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG

## VARIABLES
initial_height = 1.7 #initial height of projectile above DEM when shot [m]
alpha_min = -90.0 #minimum of vertical angle range at which the projectile is shot [°]
alpha_max = 15.0 #maximum of vertical angle range at which the projectile is shot [°]
gravitational_acceleration = -9.81 #gravitational acceleration [m/s^2]
initial_velocity = 50 #initial velocity of projectile when shot [m/s]
air_density = 1.225 #air density [kg/m^3]
drag_coefficient = 0.47 #aerodynamic drag coefficient of projectile
cross_sectional_area = 0.001963 #cross-sectional area of the projectile [m^2]
mass = 0.100 #projectile mass [kg]
dalpha = 5 #step in vertical angle range [°]
trajectory_segment_width = None #distance step, at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step
eyes_height = 1.6 #shooter eye height above DEM for viewshed [m]
target_height = 1.7 #target height for viewshed [m]
wall_height = 4.0 #obstacle/wall height (if obstacle option is used) [m]
constant = 1 #constant multipling the drag coefficient within wobble distance of an arrow
area_addition = 0.0 #average addition to cross-sectional area of an arrow within wobble distance [m^2]
wobble_distance = 40 #wobble distance - distance at which an arrow stops wobbling [m]



start = timer()


main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed, use_lines, EPSG,
         cumulative_throwshed, initial_height, initial_velocity, drag_coefficient, cross_sectional_area, mass,
         eyes_height, target_height, wall_height, constant, area_addition, wobble_distance, band_number=band_number,
         interpolation=interpolation, alpha_min=alpha_min, alpha_max=alpha_max,
         gravitational_acceleration=gravitational_acceleration, air_density=air_density, dalpha=dalpha,
         trajectory_segment_width=None)

end = timer()
print('Duration:', end - start)
