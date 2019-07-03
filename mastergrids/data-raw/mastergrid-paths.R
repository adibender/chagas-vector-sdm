## define mastergrid_paths
path_mastergrids <- "/media/z/mastergrids/"
# modis data
path_modis <- paste0(path_mastergrids, "MODIS_Global/")
# landcover
path_landcover <- paste0(path_modis, "MCD12Q1_Annual_Landcover/")
# modis reflectance
path_reflectance <- paste0(path_modis, "MCD43D6_v6_BRDF_Reflectance/")
path_reflectance_evi   <- paste0(path_reflectance, "EVI_v6/")
path_reflectance_tcb   <- paste0(path_reflectance, "TCB_v6/")
path_reflectance_tcw   <- paste0(path_reflectance, "TCW_v6/")
# modis LST
path_lst <- paste0(path_modis, "MOD11A2_v6_LST/")
path_lst_day   <- paste0(path_lst, "LST_Day/")
path_lst_night <- paste0(path_lst, "LST_Night/")
path_lst_diff  <- paste0(path_lst, "LST_DiurnalDifference/")

## Other global variables
path_other_vars <- paste0(path_mastergrids, "Other_Global_Covariates/")
# Accessibility
path_accessibility_weiss <- paste0(path_other_vars, "Accessibility/Weiss/")
# Elevation
path_elevation <- paste0(path_other_vars, "Elevation/")
# NightTimeLights
path_nighttimelights <- paste0(path_other_vars, "NightTimeLights/")
# Population
path_population <- paste0(path_other_vars, "Population/Worldpop_GPWv4_Hybrid_201708/")
# Rainfall
path_rainfall <- paste0(path_other_vars, "Rainfall/CHIRPS/")
# Urban areas
path_urban <- paste0(path_other_vars, "UrbanAreas/Global_Urban_Footprint/From_12m/1km/")

mastergrid_paths <- list(
  path_mastergrids, path_accessibility_weiss, path_elevation, path_landcover,
  path_modis, path_lst_day, path_lst_night, path_lst_diff,
  path_nighttimelights, path_other_vars, path_population,
  path_reflectance_evi, path_reflectance_tcw, path_reflectance_tcb,
  path_rainfall, path_urban)

path_names <- c(
  "mastergrid", "accessibility", "elevation", "landcover",
  "modis", "lst_day", "lst_night", "lst_diff",
  "nighttimelights", "other", "population",
  "reflectance_evi", "reflectance_tcw", "reflectance_tcb",
  "rainfall", "urban")
names(mastergrid_paths) <- path_names
use_data(mastergrid_paths, overwrite = TRUE)


landcover_df <- data.frame(
  class = paste0("Class ", c(paste0(rep(0, 10), 0:9), 10:18)),
  name = c("Water", "Evergreen Needleleaf Forest", "Evergreen Broadleaf Forest",
   "Deciduous Needleleaf Forest", "Deciduous Broadleaf Forest", "Mixed Forest",
   "Closed Shrublands", "Open Shrublands", "Woody Savannas", "Savannas",
   "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built up",
   "Cropland natural vegetation mosaic", "Snow and Ice",
   "Barren or sparsely populated", "Unclassified", "No Data"))

use_data(landcover_df, overwrite = TRUE)

landcover_queries <- purrr::map(
  purrr::map_chr(
    paste0("Class", c(paste0(rep(0, 10), 0:9), 10:18)),
    ~paste0(.x, ".*Annual.Data.5km.percentage")),
  ~ mastergrids::search_files(.x, path = mastergrid_paths[["landcover"]]))

use_data(landcover_queries, overwrite = TRUE)
