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
         interpolation=1, alpha_min=-90, alpha_max=90, gravitational_acceleration=-9.81, air_density=1.225, dalpha=5,
         trajectory_segment_width=None):
    """Just main function with controls, global variables settings and triggers to other functions"""
    # making sure the vertical angle has correct range
    if alpha_max < alpha_min:
        print("Minimal vertical shooting angle higher than maximal.")
        exit()
    # Global variables
    global SRS, DP, PLP, LLP, TOF, TF, UV, CV, BN, INT, BF, VF, TSW, AL, DB, DA, DGT, DMINH, DMAXH, IH, IV, DC, CSA, M,\
        CONST, AA, WD, GA, AD, EH, TH
    # CRS and other variable definition
    SRS = osr.SpatialReference()
    SRS.ImportFromEPSG(EPSG)
    TOF, TF, UV, CV, BN, INT, VF, IH, IV, DC, CSA, M, CONST, AA, WD, GA, AD, EH, TH = throwshed_output_folder,\
        throwshed_file, use_viewshed, cumulative_throwshed, band_number, interpolation, 'viewshed', initial_height, \
        initial_velocity, drag_coefficient, cross_sectional_area, mass, constant, area_addition, wobble_distance, \
        gravitational_acceleration, air_density, eyes_height, target_height
    # get DEM data and assign them as global
    DB, DA, DGT = get_raster_from_file(dem_path)
    # assign trajectory segment width
    TSW = np.min(np.abs([DGT[1],DGT[5]])) if trajectory_segment_width == None else trajectory_segment_width
    # obtain list of point geometries
    point_geom_list = get_geom_list_from_file(point_layer_path)
    # burn lines as obstacles into DEM, creating new DEM (Digital terrain model -> Digital surface model)
    if line_layer_path:
        burn_obstacles(line_layer_path, wall_height)
    # get minimum and maximum DEM height
    DMINH, DMAXH = get_min_max_height()
    # obtain list of vertical angles
    AL = np.arange(np.radians(alpha_min), np.radians(alpha_max + dalpha), np.radians(dalpha))
    # cycle calculating throwshed for each point
    for point_geom in point_geom_list:
        # compute throwshed for 1 point
        throwshed(point_geom)



def get_raster_from_file(file_path):
    """Get DEM array and geotransformation data"""
    # import DEM datasource
    dem_ds = gdal.Open(file_path)
    # select band
    dem_band = dem_ds.GetRasterBand(BN)
    # DEM cell values into array
    dem_array = dem_band.ReadAsArray()
    # transformation data describing DEM
    dem_gt = dem_ds.GetGeoTransform()
    dem_ds = None
    return dem_array, dem_band, dem_gt

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
    # list of geometries
    geom_list = [layer.GetFeature(number).GetGeometryRef() for number in range(0, layer.GetFeatureCount())]
    layer_ds = None
    return geom_list

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
    buffer_ds, buffer_band = create_raster_file(BFT, 0)
    # Buffer polygon is rasterized
    gdal.RasterizeLayer(buffer_ds, [1], buffer_outlayer, burn_values=[wall_height])
    # Sum of initial dem and buffer rasters
    buffer_array = buffer_band.ReadAsArray()
    global DA, DB
    DA = np.add(DA, buffer_array)
    buffer_ds = buffer_outlayer = None
    # create new raster datasource for DEM (DSM)
    dem_ds, DB = create_raster_file(BFT, 1)
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

def create_raster_file(raster_name, method):
    """Creates raster file. Method 0 returns empty datasource and band. Method 1 does the same and also writes array"""
    # create driver and output data source
    outds = gdal.GetDriverByName('GTiff').Create(TOF + "\\" + raster_name + ".tif", xsize=DA.shape[1],
                                                 ysize=DA.shape[0], bands=1, eType=gdal.GDT_Float32)
    # assign geotransformation, projection, band and nodata settings
    outds.SetGeoTransform(DGT)
    outds.SetProjection(SRS.ExportToWkt())
    raster_band = outds.GetRasterBand(1)
    raster_band.SetNoDataValue(-9999)
    if method == 1:
        raster_band.WriteArray(DA)
    return outds, raster_band

def remove_temp_files(temp_file):
    """Deletes all temporary files with assigned name"""
    for format in [file.split('.')[1] for file in os.listdir(TOF) if file.split('.')[0] == temp_file]:
        os.remove(TOF + '\\' + temp_file + "." + format)

def throwshed(point_geom):
    """Calculates throwshed for 1 point"""
    # get X and Y of shooter's place
    X_coor_point = point_geom.GetX()
    Y_coor_point = point_geom.GetY()
    # get interpolated height of point from DEM
    global point_height
    point_height = int_function(X_coor_point, Y_coor_point)
    # generate set of trajectories for vertical angle range with basic step
    trajectory_base_set()
    # CONTINUE HERE

def int_function(X, Y):
    """Interpolates height of point from DEM cells"""
    # nearest neighbour
    if INT == 0:
        column = round(np.abs((X_coor_point - (DGT[0] + DGT[1] / 2)) / DGT[1]))
        row = round(np.abs((Y_coor_point - (DGT[3] + DGT[5] / 2)) / DGT[5]))
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

def trajectory_base_set():
    """Generate set of trajectories for vertical angle range with basic step"""
    # trajectory dictionary, that will contain all trajectories, their initial shooting angle etc.
    global TS
    # dictionary's key is vertical angle value and the value (first one, in brackets) is list of heights in one trajectory
    TS = {alpha: [generate_trajectory(alpha)] for alpha in AL}

def generate_trajectory(alpha):
    """Generates trajectory from input parameters"""
    # initial drag
    d = -AD * IV ** 2 * DC * CONST * (CSA + AA) / (2 * M)
    # radial distance step
    r = TSW
    # initial projectile velocity
    V = IV
    # first/shooter point coordinates
    x = [0.0]
    y = [IH + point_height]
    # list of all trajectory points' heights
    height_list = []
    # velocity x and y elements
    V_x = V * np.cos(alpha)
    V_y = V * np.sin(alpha)
    # drag x and y elements
    d_x = d * np.cos(alpha)
    d_y = d * np.sin(alpha)
    # time step is set to value approximate to half of cell size
    dt = TSW / V / 2
    # x and y steps
    dX = V_x * dt + 1 / 2 * d_x * dt ** 2
    dY = V_y * dt + 1 / 2 * (d_y + GA) * dt ** 2
    # cycle calculating every new trajectory point one by one
    while True:
        # coords (only 2 last are kept)
        x.append(x[-1] + dX)
        y.append(y[-1] + dY)
        # when the height is finally less than the minimal DEM height, cycle breaks
        if y[-1] < DMINH:
            break
        # projectile height is interpolated every TSW-fold distance
        if x[-1] > r:
            height_list.append(((x[-1] - r) * (y[-2] - y[-1])) / (x[-1] - x[-2]) + y[-1])
            r += TSW
        # new vertical angle
        alpha = np.arctan(dY / dX)
        # new velocity
        V = ((dX / dt) ** 2 + (dY / dt) ** 2) ** (1 / 2)
        # new drag
        if (x[-1] ** 2 + (y[-1] - point_height) ** 2) ** (1 / 2) < WD:
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
        dX = V_x * dt + 1 / 2 * d_x * dt ** 2
        dY = V_y * dt + 1 / 2 * (d_y + GA) * dt ** 2
        # only two last values are utilized
        del x[0], y[0]
    return height_list

#######################################################################
## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr.tif" #path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\points1.shp"   #path to point layer
line_layer_path = None#r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp" #path to line layer, if obstacles are not to be used, set to None
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed"  #path to folder, where the file will be saved
throwshed_file = r"throwshedtest1"   #name of output throwshed file

## SETTINGS
use_viewshed = 1 #utilization of viewshed, that will clip throwshed, No = 0, Yes = 1
band_number = 1 #selected band from DEM, default = 1
interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places)
EPSG = 8353 #EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG

## VARIABLES
initial_height = 1.7 #initial height of projectile above DEM when shot [m]
alpha_min = 0.0 #minimum of vertical angle range at which the projectile is shot [°]
alpha_max = 45.0 #maximum of vertical angle range at which the projectile is shot [°]
gravitational_acceleration = -9.81 #gravitational acceleration [m/s^2]
initial_velocity = 67 #initial velocity of projectile when shot [m/s]
air_density = 1.225 #air density [kg/m^3]
drag_coefficient = 2.0 #aerodynamic drag coefficient of projectile
cross_sectional_area = 0.000050 #cross-sectional area of the projectile [m^2]
mass = 0.035 #projectile mass [kg]
dalpha = 5.0 #step in vertical angle range [°]
trajectory_segment_width = None #distance step, at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step
eyes_height = 1.6 #shooter eye height above DEM for viewshed [m]
target_height = 1.7 #target height for viewshed [m]
wall_height = 4.0 #obstacle/wall height (if obstacle option is used) [m]
constant = 1 #constant multipling the drag coefficient within wobble distance of an arrow
area_addition = 0.0 #average addition to cross-sectional area of an arrow within wobble distance [m^2]
wobble_distance = 40 #wobble distance - distance at which an arrow stops wobbling [m]


main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, use_viewshed, EPSG,
         cumulative_throwshed, initial_height, initial_velocity, drag_coefficient, cross_sectional_area, mass,
         eyes_height, target_height, wall_height, constant, area_addition, wobble_distance, band_number=1,
         interpolation=1, alpha_min=-90, alpha_max=90, gravitational_acceleration=-9.81, air_density=1.225, dalpha=5,
         trajectory_segment_width=None)


















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
