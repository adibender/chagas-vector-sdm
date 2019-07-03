## setup
library(readr)
# library(rgdal)
library(dplyr)
library(purrr)
library(furrr)
library(forcats)
library(checkmate)
# devtools::load_all("~/Dropbox/GitHub/mastergrids")
# devtools::load_all("~/Dropbox/GitHub/tcruzi/tcruziutils")
library(mastergrids)
library(tcruziutils)

## spatial infromation
# projection and extent
extent_tcruzi <- tcruziutils::extent_tcruzi

# Administrative districts data
shp_admin_2 <- rgdal::readOGR("../polygon boundaries/admin2013_2.shp") %>%
  raster::crop(extent_tcruzi) %>%
  add_projection()

shp_admin_1 <- rgdal::readOGR("../polygon boundaries/admin2013_1.shp") %>%
  raster::crop(extent_tcruzi) %>%
  add_projection()

shp_admin <- rbind(shp_admin_1, shp_admin_2)

saveRDS(shp_admin, "../preprocessing/shp_admin.Rds")
saveRDS(shp_admin_1, "../preprocessing/shp_admin_1.Rds")
saveRDS(shp_admin_2, "../preprocessing/shp_admin_2.Rds")

## country polygons
country_iso_a3 <- c("ARG", "BHS", "BLZ", "BOL", "BRA", "CAN", "CHL", "COL",
  "CRI", "CUB", "DOM", "ECU", "FRA", "GTM", "GUY", "HND", "HTI", "JAM", "MEX",
  "NIC", "PAN", "PER", "PRI", "PRY", "SLV", "SUR", "TTO", "URY", "USA", "VEN",
  "GUF")

temp_path <- tempdir()
country_polygons <- purrr::map(country_iso_a3,
  ~raster::getData("GADM", country = ., level = 0, path = temp_path))
names(country_polygons) <- purrr::map_chr(country_polygons, ~ .x@data$NAME)
saveRDS(country_polygons, "country_polygons.Rds")
countries <- do.call(rbind, country_polygons)
countries <- raster::crop(countries, tcruziutils::extent_tcruzi)
saveRDS(countries, "../preprocessing/countries.Rds")

## paths
data("tcruzi_queries",  package = "tcruziutils")

## store raster file that with proper extent and ocean water set to NA
# this raster file will be used as reference raster in the entire project, i.e.,
# create prediction data frames, visualizations, etc.
tcruzi_grid <- raster::raster(tcruzi_queries[["elevation"]]) %>%
  raster::crop(extent_tcruzi)
raster::values(tcruzi_grid)[!is.na(raster::values(tcruzi_grid))] <- 0L
raster::writeRaster(tcruzi_grid, "tcruzi_grid.tif", format = "GTiff",
  overwrite = TRUE)

# construct
future::plan("sequential")
rasters <- future_map(tcruzi_queries,
  function(query) {
    # select latest file (prediction will be made for the last available time)
    file <- query[length(query)]
    r_i  <- raster::raster(file) %>% raster::crop(raster::extent(tcruzi_grid))
    # adjust spatial resolution (all to 5x5 grids)
    if (!all(raster::res(r_i) == raster::res(tcruzi_grid))) {
      r_i <- raster::resample(r_i, tcruzi_grid)
    }
    raster::values(r_i)[is.na(raster::values(tcruzi_grid))] <- NA
    r_i
  })

env_grids <- reduce(rasters, raster::addLayer)
names(env_grids) <- names(tcruzi_queries)
env_grids <- raster::brick(env_grids, filename = "env_grids.grd",
  bylayer = FALSE, format = "raster", overwrite = TRUE)

## Vector
# vector species presence
presence_vector <- read_csv(
    file = "../infection data/combined vector occurence.csv",
    na = c("", " ", "NR", "NA")) %>%
  mutate(
    area = 25L,
    Start_year = as.integer(substr(year, 1L, 4L)),
    End_year   = as.integer(substr(year, 6L, 9L))) %>%
  rename("reference" = "associatedReferences?") %>%
  mutate(
    Public_year = as.integer(stringr::str_extract(.data$reference, "[0-9]{4}")),
    Public_year = ifelse(Public_year > 2018L, NA, Public_year),
    End_year    = ifelse(is.na(End_year), Start_year, End_year)) %>%
  select(scientificName, Start_year, End_year, Public_year, area,
    individualCount, starts_with("decimal"), habitat, reference) %>%
  rename(
    species    = scientificName,
    Latitude   = decimalLatitude,
    Longitude  = decimalLongitude,
    n_observed = individualCount) %>%
  rename_all(~tolower(sub(" ", "_", .)))
saveRDS(presence_vector, "../preprocessing/presence_vector.Rds")

## environmental grids for each year
# used to predict PPM intensity of species occurrence and sampling
# to construct weights for polygon-level infection prevalence data
files_years_df <-
  mastergrids:::group_files_by_year(tcruzi_queries, years = 2000:2018)
# construct
# future::plan("sequential") # <- change here for parallelisation over years
## WARNING: Requieres lots of RAM (in sequentiall mode ~ 32 GB)
walk(seq_len(nrow(files_years_df)),
  function(year_ind) {
    print(files_years_df$year[year_ind])
    rasters <- map(files_years_df[year_ind, "files"][[1]],
      ~raster::raster(.x) %>% raster::crop(raster::extent(tcruzi_grid)))
    # adjust spatial resolution (all to 5x5 grids)
    rasters <- map(rasters, function(r_i) {
      if (!all(raster::res(r_i) == raster::res(tcruzi_grid))) {
        r_i <- raster::resample(r_i, tcruzi_grid)
      }
      raster::values(r_i)[is.na(raster::values(tcruzi_grid))] <- NA
      r_i
    })
  env_grids_year_i <- reduce(rasters, raster::addLayer)
  names(env_grids_year_i) <- names(tcruzi_queries)
  raster::brick(
    env_grids_year_i,
    filename  = paste0("env_grids_", files_years_df$year[year_ind], ".grd"),
    bylayer   = FALSE,
    format    = "raster",
    overwrite = TRUE)
})
