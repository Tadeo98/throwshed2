#######################################################################
## THROWSHED ##
#######################################################################
# Parameters and computation start #
#######################################################################

"""
Parameters definition and calling of main throwshed function to start the computation.
"""

import throwshed2

def main():
    #######################################################################
    ## PATHS
    dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\case_study\Yann_Waersegers\project\data\dtm1m_edited.tif" #path to DEM
    point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\case_study\Yann_Waersegers\project\data\points.shp"   #path to point layer, can't be situated on raster bordering column or row (due to linear interpolation or the multiprocessing azimuth computatuions)
    line_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\case_study\Yann_Waersegers\project\data\lines.shp" #path to line layer
    throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\data\throwshed\Siege_of_Luxembourg"  #path to folder, where the file will be saved
    throwshed_file = r"defense_cannon_18pf"   #name of output throwshed file

    ## SETTINGS
    throwshed_mode = 0 #what type of throwshed will be calculated, simple safety zone (cells within safety field) = 0, regular throwshed with trajectory assessment = 1
    use_viewshed = 0 #utilization of viewshed, that will clip throwshed, No viewshed = 0, Clip by visible areas = 1, Clip by invisible areas = -1
    use_lines = 0 #utilization of line layer, where lines serve as obstacles or walls and will be burnt into DEM, No = 0, Yes = 1
    band_number = 1 #selected band from DEM, default = 1
    interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
    cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0 (binary), Yes = 1 (Apropriate with more than 1 shooting places)
    EPSG = 2169 #EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG
    atmosphere_type = 1 #standard atmospheres - Army Standard Metro (0) and ICAO (1)
    numerical_method = "euler2D" #numerical method that calculates the ballistic trajectory, string that can be euler2D (in vertical plane only), euler3D, heun2D, heun3D
    trajectory_segment_dimension = 1 #decides whether trajectory segment size stands for its width (horizontal) = 0, or length (slant) = 1 (if set to width, shooting angles close to -90 and 90 are dangerous to use as minimum and maximum of the range)
    irregular_projectile = 0 #use diameter to calculate symmetric circle cross-sectional area of projectile = 0, use cross_sectional_area for irregular shape of the projectile = 1
    cores = 12 #number of cores allocated for multiprocessing, 0-n, where 0 means all cores

    ## VARIABLES
    initial_height = 1.0 #initial height of projectile above DEM when shot [m]
    alpha_min = 0.5 #minimum of vertical angle range at which the projectile is shot [°]
    alpha_max = 11.55 #maximum of vertical angle range at which the projectile is shot [°]
    gravitational_acceleration = 9.81 #gravitational acceleration [m/s^2]
    initial_velocity = 574.82 #initial velocity of projectile when shot [m/s]
    temperature = 15 #air temperature at shooting site [°C]
    air_density = 1.225 #air density [kg/m^3]
    drag_to_mach = "Miller and Bailey 130.5mm sphere" #aerodynamic drag coefficient of projectile (constant, list or drag table)
    diameter = 0.1305 #diameter of the projectile [m]
    mass = 8.1608 #projectile mass [kg]
    cross_sectional_area = 0.0000000 #cross sectional area of the projectile, used instead of diameter when the projectile is irregular [m^2]
    trajectory_segment_size = None #distance step (length or width), at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step (equal or less than raster resolution)
    azimuth_min = 110 #minimum of the azimuth range [°]
    azimuth_max = 160 #maximum of the azimuth range [°]
    eyes_height = 1.6 #shooter eye height above DEM for viewshed [m]
    target_height = 1.7 #target height for viewshed [m]
    wall_height = 6 #obstacle/wall height or ditch depth (if obstacle option is used) [m]
    wall_width = 3 #obstacle/wall width (if obstacle option is used) [m]
    oscillation_frequency = 0 #frequency of area and drag coefficient change for projectile oscillation/rotation, 0 means no oscillation/rotation effect [Hz] or [s^-1]
    oscillation_distance = 20 #distance (of curved ballistic trajectory) up to which the rotation/oscillation of projectile linearly decreases (no rotation/oscillation afterwards), if 0, the rotation/oscillation effect remains same for whole trajectory [m]
    peak_drag = 3.56 #drag coefficient at peak of projectile's rotation/oscillation right after shooting
    peak_area = 0.0005972 #cross-sectional area at peak of projectile's rotation/oscillation right after shooting [m^2]


    throwshed2.main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed,
             use_lines, cumulative_throwshed, EPSG, atmosphere_type, numerical_method, trajectory_segment_dimension, irregular_projectile, cores,
             initial_height, initial_velocity, drag_to_mach, temperature, diameter, mass, cross_sectional_area, trajectory_segment_size, azimuth_min, azimuth_max,
             eyes_height, target_height, wall_height, wall_width, peak_drag, peak_area, oscillation_distance, oscillation_frequency, band_number=band_number,
             interpolation=interpolation, alpha_min=alpha_min, alpha_max=alpha_max,
             gravitational_acceleration=gravitational_acceleration, air_density=air_density)

if __name__ == "__main__":
    main()