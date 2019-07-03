data("mastergrid_paths", package = "mastergrids")

tcruzi_paths <- mastergrid_paths[c(
  "accessibility", "elevation", "landcover",
  "lst_day", "lst_night", "lst_diff",
  "nighttimelights", "population",
  "reflectance_evi", "reflectance_tcb", "reflectance_tcw",
  "rainfall", "urban")]

tcruzi_patterns <- paste0(c(
  "accessibility_to_cities_2015.*", ".*elevation.*5km.mean",
  ".*Landcover_Class02.*Annual.Data.5km.percentage",
  rep("LST.*Annual.mean.5km.mean", 3),
  "VIIRS.*Annual.5km.MEAN*", ".*Pop.*5km.*",
  "EVI.*Annual.mean.*5km.mean", "TCB.*Annual.mean.*5km.mean",
  "TCW.*Annual.mean.*5km.mean",
   "chirps-v2.*Annual.sum.5km.*", "Global.*Unclipped"), ".tif$")

usethis::use_data(tcruzi_paths, tcruzi_patterns, overwrite = TRUE)


tcruzi_queries <- purrr::map2(
  tcruzi_paths[c(1, rep(2, 2), rep(3, 19), 4:13)],
  c("accessibility_to_cities_2015.*", ".*elevation.*5km.mean",
  "*SlopePCT_Corrected",
  paste0(paste0("Class", c(paste0(rep(0, 10), 0:9), 10:18)),
    ".*Annual.Data.5km.percentage"),
  rep("LST.*Annual.mean.5km.mean", 3),
  "VIIRS.*Annual.5km.MEAN*", ".*Pop.*5km.*",
  "EVI.*Annual.mean.*5km.mean", "TCB.*Annual.mean.*5km.mean",
  "TCW.*Annual.mean.*5km.mean",
   "chirps-v2.*Annual.sum.5km.*", "Global.*Unclipped"),
  ~ mastergrids::search_files(.y, path = .x))

names(tcruzi_queries) <- c("accessibility", "elevation", "slope",
 paste0("landcover", c(paste0(rep(0, 10), 0:9), 10:18)),
  "lst_day", "lst_night", "lst_diff",
  "nighttimelights", "population",
  "reflectance_evi", "reflectance_tcb", "reflectance_tcw",
  "rainfall", "urban")

usethis::use_data(tcruzi_queries, overwrite = TRUE)
