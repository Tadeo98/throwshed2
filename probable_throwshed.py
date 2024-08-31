#######################################################################
## THROWSHED ##
#######################################################################
# Probable throwshed with Monte Carlo simulation #
#######################################################################

"""
Probable throwshed with Monte Carlo simulation based on deviations/range of input physical parameter or DEM.
Parameters are arrays of random normally distributed or within range uniformly distributed values around the mean.
Uniformly distributed values within the range is possible only for physical parameters, DEM errors will be distributed
normally only.
Parameter drag_to_mach can also hold non-zero deviation but only if it is meant to be a constant, because it can bear
list values anyway (specific drag values for specific Mach, often used for supersonic projectiles).
DEM raster can be recalculated by adding rasters with deviations/errors to it to serve as another mean of the simulation.
Computation of individual throwsheds is conducted many times and these are all added into one probable throwshed file.
Option with multiple shooters is possible but only with the cumulative mode being off (setting unreachable).
If the parameter is considered deterministic - same for all computations, set the deviation to 0.0.
"""

import throwshed2
import os
import numpy as np
from osgeo import gdal, osr

#######################################################################

#FUNCTIONS

def main():
    global SRS, dds, db, da, dgt, dndv
    throwshed_file = r"throwshed_temp"  # name of output throwshed file
    cumulative_throwshed = 0  # Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places) --- For Monte Carlo simulation this has to be set to 0
    # name of the temporary error DEM file
    error_DEM = "error_DEM"
    # path to DEM is changed to the one with errors
    DEM_path = throwshed_output_folder + "\\" + error_DEM + ".tif"
    # a short message that for DEM only normal distribution is possible, calculation is conducted anyway
    if MCS_dist == 0:
        print("Distribution is set to uniform, but for DEM, normal distribution will be used.")
    # program exits when the drag is set as a list and is chosen to be generated randomly at the same time
    if drag_to_mach[1] and (type(drag_to_mach[0]) == str or type(drag_to_mach[0]) == list):
        print("Drag cannot be chosen as a randomly generated array because it is already set as a list of drag to mach function. This works only with constant drag.")
        exit()
    # CRS definition
    SRS = osr.SpatialReference()
    SRS.ImportFromEPSG(EPSG)
    # get DEM parameters
    dds, db, da, dgt, dndv = get_raster_from_file(dem_path,[1])
    # create raster arrays with zeroes, Probable Throwshed Raster (Array)
    PTRA = [np.zeros((da[0].shape[0], da[0].shape[1]), np.int16), np.zeros((da[0].shape[0], da[0].shape[1]), np.int16)]
    # put the parameters (all to which an array can be assigned) into dictionary
    parameters = {
        "initial_height": initial_height,
        "gravitational_acceleration": gravitational_acceleration,
        "initial_velocity": initial_velocity,
        "temperature": temperature,
        "air_density": air_density,
        "diameter": diameter,
        "cross_sectional_area": cross_sectional_area,
        "mass": mass,
        "drag_to_mach": drag_to_mach
    }
    # fill the parameters with randomly generated normally distributed or uniformly distributed (with a constant step) values in arrays
    if MCS_dist:
        for key, value in parameters.items():
            parameters[key] = np.random.normal(value[0],value[1],MCS_number)
    else:
        for key, value in parameters.items():
            parameters[key] = np.linspace(value[0]-value[1],value[0]+value[1],MCS_number)
    # cycle that creates DEM with errors, generates throwsheds with all the different randomly generated values of parameters
    for i in range(MCS_number):
        # create Error Raster (Array) with set DEM error
        ERA = np.random.normal(0, DEM_dev, (da[0].shape[0], da[0].shape[1]))
        # create mask with False at positions where NoData values lie
        mask = (da[0] != dndv)
        # New DEM Array by adding ERA to original DEM array, NoData values unchanged
        NDEMA = [np.where(mask, da[0] + ERA, dndv)]
        # create temporary error DEM raster file
        create_raster_file(error_DEM, NDEMA, gdal.GDT_Float32, dndv, [1])
        # generate one throwshed
        throwshed2.main(DEM_path,point_layer_path,line_layer_path,throwshed_output_folder,throwshed_file,
        throwshed_mode,use_viewshed,use_lines,cumulative_throwshed,EPSG,atmosphere_type,numerical_method,
        trajectory_segment_dimension,irregular_projectile,cores,parameters.get("initial_height")[i],parameters.get("initial_velocity")[i],
        parameters.get("drag_to_mach")[i],parameters.get("temperature")[i],parameters.get("diameter")[i],
        parameters.get("mass")[i],parameters.get("cross_sectional_area")[i],azimuth_min,azimuth_max,target_height,wall_height,
        wall_width,peak_drag,peak_area,oscillation_distance,oscillation_frequency,band_number=band_number,
        interpolation=interpolation,alpha_min=alpha_min,alpha_max=alpha_max,
        gravitational_acceleration=parameters.get("gravitational_acceleration")[i],
        air_density=parameters.get("air_density")[i],trajectory_segment_size=None)
        # get throwshed parameters
        tds, tb, ta, tgt, tndv = get_raster_from_file(os.path.join(throwshed_output_folder, throwshed_file + ".tif"),[1,2])
        # create mask with False at positions where NoData values lie
        mask1, mask2 = (ta[0] != dndv), (ta[1] != dndv)
        # add new throwshed raster, NoData values unchanged
        PTRA = [np.where(mask1, PTRA[0] + ta[0], dndv), np.where(mask2, PTRA[1] + ta[1], dndv)]
        tds = tb = ta = tgt = tndv = None
        # remove temporary error DEM file(s)
        throwshed2.remove_temp_files(error_DEM)
        # just a note
        print(f"Throwshed number {i+1} was added.")
    # create mask with False at positions where NoData values lie
    mask1, mask2 = (PTRA[0] != dndv), (PTRA[1] != dndv)
    # divide values by number of repetitions and get probability values ranging from 0 to 1, NoData values unchanged
    PTRA = [np.where(mask1, PTRA[0]/MCS_number, dndv), np.where(mask2, PTRA[1]/MCS_number, dndv)]
    # create probable throwshed raster file
    create_raster_file(PTF, PTRA, gdal.GDT_Float32, dndv, [1,2])
    # close all DEM parameters
    dds = db = da = dgt = dndv = None
    # remove temporary throwshed file
    os.remove(throwshed_output_folder + '\\' + throwshed_file + ".tif")

def get_raster_from_file(file_path, bands):
    """Get DEM datasource, band, array and geotransformation data and nodata value"""
    # import DEM datasource
    dem_ds = gdal.Open(file_path)
    dem_band, dem_array = [], []
    for band in bands:
        # select band
        dem_band.append(dem_ds.GetRasterBand(band))
        # DEM cell values into array
        dem_array.append(dem_band[band-1].ReadAsArray())
    # transformation data describing DEM
    dem_gt = dem_ds.GetGeoTransform()
    # nodata value
    no_data_value = dem_band[0].GetNoDataValue()
    return dem_ds, dem_band, dem_array, dem_gt, no_data_value

def create_raster_file(raster_name, dem_array, GDT, no_data, bands):
    """Creates raster file out of arrays."""
    # create driver and output data source
    outds = gdal.GetDriverByName('GTiff').Create(throwshed_output_folder + "\\" + raster_name + ".tif", xsize=da[0].shape[1],
                                                 ysize=da[0].shape[0], bands=len(bands), eType=GDT)
    # assign geotransformation, projection, band and nodata settings
    outds.SetGeoTransform(dgt)
    outds.SetProjection(SRS.ExportToWkt())
    for band in bands:
        raster_band = outds.GetRasterBand(band)
        raster_band.SetNoDataValue(no_data)
        raster_band.WriteArray(dem_array[band-1])

#######################################################################

## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr_clip.tif"  # path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\point.shp"  # path to point layer
line_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp"  # path to line layer
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed\arrow_wobble\AVIb0"  # path to folder, where the file will be saved
PTF = r"probable_velocity_normal_n20_s5_single_simple_f60_d20"  # name of output Probable Throwshed File

## SETTINGS
throwshed_mode = 0  # what type of throwshed will be calculated, simple safety zone (cells within safety field) = 0, regular throwshed with trajectory assessment = 1
use_viewshed = 0 #utilization of viewshed, that will clip throwshed, No viewshed = 0, Clip by visible areas = 1, Clip by invisible areas = -1
use_lines = 0  # utilization of line layer, where lines serve as obstacles or walls and will be burnt into DEM, No = 0, Yes = 1
band_number = 1  # selected band from DEM, default = 1
interpolation = 0  # interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
EPSG = 8353  # EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG
atmosphere_type = 1  # standard atmospheres - Army Standard Metro (0) and ICAO (1)
numerical_method = "euler2D"  # numerical method that calculates the ballistic trajectory, string that can be euler2D (in vertical plane only), euler3D, heun2D, heun3D
trajectory_segment_dimension = 1  # decides whether trajectory segment size stands for its width (horizontal) = 0, or length (slant) = 1 (if set to width, shooting angles close to -90 and 90 are dangerous to use as minimum and maximum of the range)
irregular_projectile = 1 #use diameter to calculate symmetric circle cross-sectional area of projectile = 0, use cross_sectional_area for irregular shape of the projectile = 1
MCS_dist = 1 # decides whether values generate as uniformly distributed within the range (MCS_dev_range) = 0, or as normally distributed due to parameter deviation (MCS_dev_range) = 1
cores = 0 #number of cores allocated for multiprocessing, 0-n, where 0 means all cores

## VARIABLES
#following 8 can be set as an array of randomly generated values, format [mean, deviation/one way range]:
initial_height = [1.7, 0.0] #initial height of projectile above DEM when shot [m]
gravitational_acceleration = [9.810, 0.000] #gravitational acceleration [m/s^2]
initial_velocity = [50.0, 5.0] #initial velocity of projectile when shot [m/s]
temperature = [15.0, 0.0] #air temperature at shooting site [°C]
air_density = [1.225, 0.000] #air density [kg/m^3]
drag_to_mach = [0.47, 0.00] #aerodynamic drag coefficient of projectile (constant, list or drag table)
diameter = [0.050, 0.000] #diameter of the projectile [m^2]
cross_sectional_area = [0.0000938, 0.0000000] #cross sectional area of the projectile, used instead of diameter when the projectile is irregular [m^2]
mass = [0.100, 0.000] #projectile mass [kg]
#other variables:
alpha_min = -90.0 #minimum of vertical angle range at which the projectile is shot [°]
alpha_max = 90.0 #maximum of vertical angle range at which the projectile is shot [°]
trajectory_segment_size = None  # distance step (length or width), at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step (equal or less than raster resolution)
azimuth_min = 45 #minimum of the azimuth range [°]
azimuth_max = 180 #maximum of the azimuth range [°]
eyes_height = 1.6  # shooter eye height above DEM for viewshed [m]
target_height = 1.7  # target height for viewshed [m]
wall_height = 4.0  # obstacle/wall height or ditch depth (if obstacle option is used) [m]
wall_width = 0.2  # obstacle/wall width (if obstacle option is used) [m]
oscillation_frequency = 60 #frequency of area and drag coefficient change for projectile oscillation/rotation, 0 means no oscillation/rotation effect [Hz] or [s^-1]
oscillation_distance = 20 #distance (of curved ballistic trajectory) up to which the rotation/oscillation of projectile linearly decreases (no rotation/oscillation afterwards), if 0, the rotation/oscillation effect remains same for whole trajectory [m]
peak_drag = 3.56 #drag coefficient at peak of projectile's rotation/oscillation right after shooting
peak_area = 0.0005972 #cross-sectional area at peak of projectile's rotation/oscillation right after shooting [m^2]
# Other probable throwshed/Monte Carlo specifications:
MCS_number = 20 # number of Monte Carlo simulations
DEM_dev = 0.00 # DEM elevation deviation for the Monte Carlo simulation

main()