##","boundary bos for tcruzi analysis
extent_tcruzi <- raster::extent(-125, -34, -42, 46.2) # see "endemic zone" folder
usethis::use_data(extent_tcruzi, overwrite = TRUE)
bbox_tcruzi <- tmaptools::bb(extent_tcruzi,
  current.projection = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
usethis::use_data(bbox_tcruzi, overwrite = TRUE)


## base plot endemic zone
data("World", package = "tmap")

endemic_zone <- sf::st_transform(World,
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  tmaptools::crop_shape(bbox_tcruzi)
p_endemic_zone <- tmap::tm_shape(endemic_zone) + tmap::tm_borders()
usethis::use_data(p_endemic_zone, overwrite = TRUE)

## standard projection
proj_tcruzi <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
usethis::use_data(proj_tcruzi, overwrite = TRUE)

# usethis::use_data(country_polygons, overwrite = TRUE)
