#######################################################################
## THROWSHED ##
#######################################################################
# Parameters and computation start #
#######################################################################

"""
Parameters definition and calling of main throwshed function to start the computation.
"""

import throwshed2

#######################################################################
## PATHS
dem_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\dem\dmr_clip.tif" #path to DEM
point_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\point\point.shp"   #path to point layer
line_layer_path = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\line\lines1.shp" #path to line layer
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term11\Dizertacna_praca\Throwshed2\data\throwshed"  #path to folder, where the file will be saved
throwshed_file = r"throwshedtest"   #name of output throwshed file

## SETTINGS
throwshed_mode = 0 #what type of throwshed will be calculated, simple safety zone (cells within safety field) = 0, regular throwshed with trajectory assessment = 1
use_viewshed = 0 #utilization of viewshed, that will clip throwshed, No = 0, Yes = 1
use_lines = 0 #utilization of line layer, where lines serve as obstacles or walls and will be burnt into DEM, No = 0, Yes = 1
band_number = 1 #selected band from DEM, default = 1
interpolation = 0 #interpolation of DEM to calculate altitude of shooting point or compare points within the DEM-to-trajectory comparison function, Nearest neighbour = 0, Bilinear = 1
cumulative_throwshed = 1 #Calculate cumulative throwshed? No = 0, Yes = 1 (Apropriate with more than 1 shooting places)
EPSG = 8353 #EPSG code for CRS of output throwshed layer and other temporary results, must be same as DEM's EPSG
atmosphere_type = 1 #standard atmospheres - Army Standard Metro (0) and ICAO (1)
numerical_method = "euler2D" #numerical method that calculates the ballistic trajectory, string that can be euler2D (in vertical plane only), euler3D, heun2D, heun3D
trajectory_segment_dimension = 1 #decides whether trajectory segment size stands for its width (horizontal) = 0, or length (slant) = 1 (if set to width, shooting angles close to -90 and 90 are dangerous to use as minimum and maximum of the range)

## VARIABLES
initial_height = 1.7 #initial height of projectile above DEM when shot [m]
alpha_min = -90.0 #minimum of vertical angle range at which the projectile is shot [°]
alpha_max = 0.0 #maximum of vertical angle range at which the projectile is shot [°]
gravitational_acceleration = 9.81 #gravitational acceleration [m/s^2]
initial_velocity = 200 #initial velocity of projectile when shot [m/s]
temperature = 15 #air temperature at shooting site [°C]
air_density = 1.225 #air density [kg/m^3]
drag_to_mach = 0.47 #aerodynamic drag coefficient of projectile (constant, list or drag table)
diameter = 0.05 #diameter of the projectile [m^2]
mass = 0.100 #projectile mass [kg]
trajectory_segment_size = None #distance step (length or width), at which trajectory's points will be saved and compared to DEM [m], None = adjusted to DEM resolution (cell's size), any float/int value = customized distance step (equal or less than raster resolution)
eyes_height = 1.6 #shooter eye height above DEM for viewshed [m]
target_height = 1.7 #target height for viewshed [m]
wall_height = 4.0 #obstacle/wall height (if obstacle option is used) [m]
wall_width = 0.2 #obstacle/wall width (if obstacle option is used) [m]
# constant = 1 #constant multipling the drag coefficient within wobble distance of an arrow
# area_addition = 0.0 #average addition to cross-sectional area of an arrow within wobble distance [m^2]
# wobble_distance = 40 #wobble distance - distance at which an arrow stops wobbling [m]


throwshed2.main(dem_path, point_layer_path, line_layer_path, throwshed_output_folder, throwshed_file, throwshed_mode, use_viewshed,
         use_lines, cumulative_throwshed, EPSG, atmosphere_type, numerical_method, trajectory_segment_dimension,
         initial_height, initial_velocity, drag_to_mach, temperature, diameter, mass,
         eyes_height, target_height, wall_height, wall_width, band_number=band_number,
         interpolation=interpolation, alpha_min=alpha_min, alpha_max=alpha_max,
         gravitational_acceleration=gravitational_acceleration, air_density=air_density, trajectory_segment_size=None)
