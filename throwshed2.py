#######################################################################
## THROWSHED ##
#######################################################################

import os
import numpy as np
from osgeo import gdal, ogr, osr


#######################################################################
## FUNCTIONS

def main(dem_path, point_layer_path, throwshed_output_folder, throwshed_file, use_viewshed, EPSG, cumulative_throwshed,
         h, V_0, C_d, A, m, h_E, h_T, wall_height, const, A_p, w_d, line_layer_path=None,
         viewshed_file='viewshed', buffer_file='buffer', band_number=1, interpolation=1, keep_crs=1, alfa_min=-90,
         alfa_max=90, g=-9.81, ro=1.225, dalfa=5, dr=None):
    """Just main function that triggers everything else"""
    dem_array, dem_band, dem_gt = get_raster_from_file(dem_path,band_number)
    dem_height = get_min_height(dem_band)
    point_layer, point_feat_count = get_layer_from_file(point_layer_path)


def get_raster_from_file(file_path,band_number):
    """Get DEM array and geotransformation data"""
    # import DEM datasource
    dem_ds = gdal.Open(file_path)
    # select band
    dem_band = dem_ds.GetRasterBand(band_number)
    # DEM cell values into array
    dem_array = dem_band.ReadAsArray()
    # transformation data describing DEM
    dem_gt = dem_ds.GetGeoTransform()
    return dem_array, dem_band, dem_gt

def get_min_height(dem_band):
    """Get minimal DEM height"""
    # sometimes the function returns None, therefore statistics need to be calculated first
    if dem_band.GetMinimum() == None:
        dem_band.ComputeStatistics(0)
    return dem_band.GetMinimum()

def get_layer_from_file(file_path):
    """Get layer from vector layer file"""
    # import vector layer datasource
    layer_ds = ogr.Open(file_path, 0)  # 1 = editing, 0 = read only. Datasource
    # vector layer
    layer = layer_ds.GetLayer()
    # get feature count
    feat_count = layer.GetFeatureCount()
    return layer, feat_count

    # Making sure the azimuth minumum and maximum are within the range of one circle (0;2pi) and their conversion to rads
    min_azimuth, max_azimuth, dazimuth = min_azimuth / 360 % 1 * (2 * np.pi), max_azimuth / 360 % 1 * (
                2 * np.pi), np.radians(dazimuth)
    # in case of min azimuth being greater or equal than max azimuth, 2pi is added to max azimuth
    if min_azimuth >= max_azimuth:
        max_azimuth += np.pi * 2
        # kontrola, ci su hodnoty uhla vystrelu spravne definovane
        if alfa_max < alfa_min:
            print("Minimalny uhol vystrelu ma vacsiu hodnotu ako maximalny. Opravit.")
            exit()

#######################################################################
## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr.tif" #path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\points1.shp"   #path to point layer
line_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp" #path to line layer, if obstacles are not to be used, set to None
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed"  #path to folder, where the file will be saved
throwshed_file = r"throwshedtest1"   #name of output throwshed file
viewshed_file = r"viewshed" #name of temporary file containing viewshed
buffer_file = r"buffer"   #name of temporary file containing buffer

## SETTINGS
use_viewshed = 1 #utilization of viewshed, that will clip throwshed, No = 0, Yes = 1
band_number = 1 #selected band from DEM, default = 1
interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
keep_crs = 0 #use CRS from input DEM for throwshed layer, Yes = 1, No = 0 (requires setting of EPSG variable)
cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places)
EPSG = 8353 #EPSG code for CRS of output throwshed layer

## VARIABLES
h = 1.7 #initial height of projectile above DEM when shot [m]
alfa_min = 0.0 #minimum of vertical angle range at which the projectile is shot [°]
alfa_max = 45.0 #maximum of vertical angle range at which the projectile is shot [°]
g = -9.81 #gravitational acceleration [m/s^2]
V_0 = 67 #initial velocity of projectile when shot [m/s]
ro = 1.225 #air density [kg/m^3]
C_d = 2.0 #aerodynamic drag coefficient of projectile
A = 0.000050 #cross-sectional area of the projectile [m^2]
m = 0.035 #projectile mass [kg]
dalfa = 5.0 #step in vertical angle range [°]
dr = None #distance step, at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step
h_E = 1.6 #shooter eye height above DEM for viewshed [m]
h_T = 1.7 #target height for viewshed [m]
wall_height = 4.0 #obstacle/wall height (if obstacle option is used) [m]
const = 1 #constant multipling the drag coefficient within wobble distance of an arrow
A_p = 0.0 #average addition to cross-sectional area of an arrow within wobble distance [m^2]
w_d = 40 #wobble distance - distance at which an arrow stops wobbling [m]


main(dem_path, point_layer_path, throwshed_output_folder, throwshed_file, use_viewshed, EPSG, cumulative_throwshed,
    h, V_0, C_d, A, m, h_E, h_T, wall_height, const, A_p, w_d, line_layer_path=None,
    viewshed_file='viewshed', buffer_file='buffer', band_number=1, int_compare=1, keep_point_crs=1, alfa_min=-90,
    alfa_max=90, g=-9.81, ro=1.225, dalfa=5, dr=None)
















#######################################################################
## VYPOCET DMP NA ZAKLADE VOLBY POUZITIA HRADIEB V RASTRI
# Vypocet na zaklade volby pouzitia hradieb v rastri
if use_line == 0:
    # nepouziju sa hradby, len sa pripocita vyska ochodze k vyske vystreleneho sipu, resp. pripocita sa vyska 0 m ak stoji obranca/utocnik priamo na terene
    h = h+feet_height
# pouzije sa vektorova vrstva s hradbami
elif use_line == 1:
    # zatvorime pasmo aby sme don mohli dat nove (ako DMP), pretoze este sa pasmo pouziva pri tvorbe viewshedu
    dem_band = None
    # import liniovej vrstvy
    line_ds = ogr.Open(line_layer_path, 0) #1 = editing, 0 = read only. Datasource
    # liniova vrstva
    line_layer = line_ds.GetLayer()
    # ziskanie poctu linii vo vrstve
    line_count = line_layer.GetFeatureCount()

    # nastavenia polygonovej vrstvy, ktora bude obsahovat buffer
    shpdriver = ogr.GetDriverByName("ESRI Shapefile")
    buffer_outds = shpdriver.CreateDataSource(throwshed_output_folder + "\\" + buffer_file + "_temp.shp")
    srs = line_layer.GetSpatialRef()
    buffer_outlayer = buffer_outds.CreateLayer(buffer_file + "_temp", srs)
    buffer_feature = ogr.Feature(buffer_outlayer.GetLayerDefn())

    buffer_dist = dr/2 + (dem_gt[1]**2+dem_gt[5]**2)**(1/2)/2 - (dr/2)/(dem_gt[1]**2+dem_gt[5]**2)**(1/2)%1 #sirka buffera z kroku dr a uhlopriecky bunky rastra s DMR
    # cyklus citajuci vsetky linie z liniovej vrstvy
    for line_number in range(0,line_count):
        # ziskanie konkretnej linie
        line_feature = line_layer.GetFeature(line_number)
        # ziskanie geometrie konktretnej linie
        line_geom = line_feature.GetGeometryRef()
        # tvorba bufferov
        if line_number == 0:
            # ak je iba jedna linia, zapise sa sama
            buffer_geom = line_geom.Buffer(buffer_dist,1) # 1 - pocet bodov tvoriacich 90° obluk okolo koncovych bodov
        else:
            # ak je viacero linii/bufferov, tak sa pridruzia spolu cez Union do jedneho prvku (aby sa buffre neprekryvali, akoby Dissolve)
            buffer_geom = buffer_geom.Union(line_geom.Buffer(buffer_dist,1))

    # zapis spojenych bufferov do vrstvy
    buffer_feature.SetGeometry(buffer_geom)
    buffer_outlayer.CreateFeature(buffer_feature)

    # nastavenia rastra s bufferom
    buffer_ds = gdal.GetDriverByName('GTiff').Create(throwshed_output_folder + "\\" + buffer_file + "_temp.tif", xsize = dem_array.shape[1], ysize = dem_array.shape[0], bands = 1, eType = gdal.GDT_Float32)
    buffer_ds.SetGeoTransform(dem_gt)
    buffer_ds.SetProjection(srs.ExportToWkt())   #SS bude nastaveny ako bol aj pri polygonovej vrstve
    buffer_band = buffer_ds.GetRasterBand(1)
    buffer_band.SetNoDataValue(0)
    # buffer polygon sa rasterizuje, [1] - priradenie hodnot do pasma 1, burn_values=[1] - priradenie hodnot buniek = 1
    gdal.RasterizeLayer(buffer_ds, [1], buffer_outlayer, burn_values=[wall_height])

    # suma dvoch rastrov
    buffer_array = buffer_band.ReadAsArray()
    dem_array = np.add(dem_array,buffer_array)

    # vytvorenie a nastavenie rastra DMP
    dmp_driver = gdal.GetDriverByName("GTiff")
    dmp_outds = dmp_driver.Create(throwshed_output_folder + "\\" + buffer_file + "_dmp_temp.tif", xsize = dem_array.shape[1], ysize = dem_array.shape[0], bands = 1, eType = gdal.GDT_Float32)
    dmp_outds.SetGeoTransform(dem_gt)
    dmp_outds.SetProjection(srs.ExportToWkt())
    dem_band = dmp_outds.GetRasterBand(1) #sice dem_band, no spravne by malo byt dmp_band
    dem_band.WriteArray(dem_array)
    dem_band.SetNoDataValue(-9999)

    # zavretie vsetkeho okrem rastra
    buffer_ds = buffer_outds = buffer_outlayer = buffer_feature = None

    # vymazanie .shp a .tif vrstvy s bufferom
    for format in [file.split('.')[1] for file in os.listdir(throwshed_output_folder) if file.split('.')[0] == buffer_file + "_temp"]:
        os.remove(throwshed_output_folder + '\\' + buffer_file + "_temp." + format)
else:
    print("Nespravne nastavena volba uvazovania s hradbami.")
    exit()

#######################################################################
## VYPOCET THROWSHEDU Z KAZDEHO BODU

point_number_once = 0
# cyklus v ktorom sa vytvori raster throwshedu pre kazdy bod bodovej vrstvy
for point_number in range(0,point_count):
    # ziskanie konkretneho bodu
    point_feature = point_layer.GetFeature(point_number)
    # ziskanie geometrie bodu, X a Y suradnice (odpovedajuce smeru a orientacii osi v QGISe)
    point_geom = point_feature.GetGeometryRef()
    X_coor_point = point_geom.GetX()
    Y_coor_point = point_geom.GetY()
    
    # zistenie vysky bunky, na ktorej sa nachadza bod, zatial len metoda najblizsieho suseda
    dem_cell_column = round(np.abs((X_coor_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1]))
    dem_cell_row = round(np.abs((Y_coor_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5]))
    dem_cell_height = dem_array[dem_cell_row][dem_cell_column]


    # VYPOCET TRAJEKTORIE PROJEKTILU (x,y,y_r)
    #Vytvorenie listu so vsetkymi hodnotami uhla vystrelu, ktore sa pouziju v cykle
    alfa_list = np.arange(np.radians(alfa_min),np.radians(alfa_max+dalfa),np.radians(dalfa))

    y_r = []    #buduca matica s vyskami projektilu kazdych dr metrov (hodnoty v riadku vedla seba) pod kazdym uhlom vystrelu (niekolko riadkov pod sebou)
    #cyklenie sa vsetkymi hodnotami alfa
    for alfa in alfa_list:
        # Pociatocny drag
        d = -ro*V_0**2*C_d*const*(A+A_p)/(2*m)    #presny vztah
        # d = -C_d/1000*V_0**2    #priblizny vztah, pre sip postacuje

        r = dr   #radialna vzdialenost (interpolacny skok)
        V = V_0 #rychlost
        x = [0.0]   #suradnica x v 2D grafe
        # y = [h] #suradnica y v 2D grafe, iba vyska nad povrchom
        y = [h+dem_cell_height] #suradnica y v 2D grafe
        y_r1 = []    #bude obsahovat vysku kazdych dr metrov (iba jeden riadok)

        #zlozky pociatocnej rychlosti v smeroch x a y
        V_x = V*np.cos(alfa)
        V_y = V*np.sin(alfa)

        #zlozky pociatocneho dragu v smeroch x a y
        d_x = d*np.cos(alfa)
        d_y = d*np.sin(alfa)

        #kroky v x a y
        dX = V_x*dt+1/2*d_x*dt**2
        dY = V_y*dt+1/2*(d_y+g)*dt**2

        while True:
            #suradnice
            x.append(x[-1] + dX)
            y.append(y[-1] + dY)

            #when the height is finally less than the minimal DEM height, last value of recorded trajectory is set to that value
            if y[-1] < min_height:
                y_r1.append(min_height)
                break

            #v kazdom nasobku nastaveneho skoku radialnej vzialenosti sa vypocita vyska sipu
            if x[-1] > r:
                y_r1.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
                r += dr

            #novy uhol
            alfa = np.arctan(dY/dX)
            
            #nova rychlost
            V = ((dX/dt)**2+(dY/dt)**2)**(1/2)
            
            #novy drag
            if (x[-1]**2+(y[-1]-dem_cell_height)**2)**(1/2) < w_d:
                d = -ro*V**2*C_d*const*(A+A_p)/(2*m)
            else:
                d = -ro*V**2*C_d*A/(2*m)    #presny vztah
            #d = -C_d/1000*V**2    #priblizny vztah, pre sip postacuje
            
            #zlozky pociatocneho dragu v smeroch x a y
            d_x = d*np.cos(alfa)
            d_y = d*np.sin(alfa)
            #zlozky rychlosti v smeroch x a y
            V_x += d_x*dt
            V_y += (d_y+g)*dt
            #kroky v x a y
            dX = V_x*dt+1/2*d_x*dt**2
            dY = V_y*dt+1/2*(d_y+g)*dt**2

            del x[0], y[0]  # only two last values are utilized

        y_r.append(y_r1)


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
