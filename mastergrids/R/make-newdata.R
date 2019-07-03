#' Creates a data frame from grid (file)
#'
#' Creates a data frame from grid file, having as many rows as unique coordinates
#' in the grid file, excluding cells with \code{NA} values.
#'
#' @inheritParams get_grid_data
#' @importFrom raster as.data.frame
#' @importFrom dplyr rename select
#' @export
grid_to_df <- function(
  grid,
  extent     = NULL,
  mask       = NULL,
  spat_form  = ~ longitude + latitude,
  na.rm = TRUE) {

  if (testFile(grid, access = "r", extension = "tif")) {
    grid <- read_grid(grid, extent, mask)
  } else {
    grid <- grid %>% raster::crop(get_extent(extent, mask))
    if(!is.null(mask)) {
      grid <- grid %>% raster::mask(mask)
    }
  }

  grid %>%
    raster::as.data.frame(xy = TRUE, na.rm = na.rm) %>%
    rename(
      !!all.vars(spat_form)[[1]] := "x",
      !!all.vars(spat_form)[[2]] := "y")

}

#' @rdname grid_to_df
grid_as_coord_df <- function(
  grid,
  extent     = NULL,
  mask       = NULL,
  spat_form  = ~ longitude + latitude) {

  grid_to_df(grid, extent, mask, spat_form)  %>%
    select(all.vars(spat_form))

}

#' @importFrom purrr reduce
#' @importFrom furrr future_map
combine_extent <- function(
  polygon_list,
  max_extent,
  .progress = interactive()) {

  poly_extent <- future__map(
    polygons[countries],
    ~raster::intersect(raster::extent(.x), max_extent),
    .progress = .progress) %>%
  reduce(~raster::extend)

}
