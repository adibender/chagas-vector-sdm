#' Adds environmental variables from grid data to point data
#'
#' @param data A data frame with spatial coordinates
#' @param grid A list of gridded data sets. Length 1 if synoptic data, otherwise
#' one element for each year for which grid data is available.
#' Must be a named list. A column with the name of the list will be
#' added to \code{data}, containing respective grid data information for each
#' year available in \code{grid}. If \code{data} contains years that are not
#' available in \code{grid}, the first element of \code{grid} will be used.
#' @param col_name The name of the newly created variable can be specified here,
#' otherwise will be extracted from \code{names(grid)}.
#' @param year_var The name of the variable in \code{data} from which the
#' years will be extracted for which matches will be looked for in \code{grid}.
#' @inheritParams spat_plot
#'
#' @import dplyr
#' @importFrom purrr map_dfr
#' @importFrom raster raster extract
#'
#' @export
add_grid_data <- function(
  data,
  grid,
  col_name = NULL,
  year_var = "start_year",
  spat_form = ~ Longitude + Latitude) {

  if (is.null(col_name)) {
    col_name <- gsub("[0-9]", "", names(grid)[1])
  }

  long_var <- all.vars(spat_form)[[1]]
  lat_var <- all.vars(spat_form)[[2]]

  years <- as.numeric(str_extract_year(names(grid)))
  years_data <- sort(unique(data[[year_var]]))
  years_ind <- years %in% years_data
  years <- intersect(years_data, years)
  if (any(years_ind)) {
    grid <- grid[years_ind]
  } else {
    grid <- grid[1]
  }

  data_with_sp_cov <- data %>%
    split(f = data[[year_var]]) %>%
    map_dfr(function(z) {
      year <- unique(z[[year_var]])
      # print(year)
      ind_grid <- which(years == year)
      if (length(ind_grid) == 0) {
        ind_grid <- 1
      }
      dat_sp <- z %>%
        select(one_of(c("id", long_var, lat_var))) %>%
        filter(!is.na(.data[[lat_var]]))
      coordinates(dat_sp) <- spat_form
      dat_sp[[col_name]] <- extract(raster(grid[[ind_grid]]),
       dat_sp)
      dat_sp %>% as.data.frame()
    })

  data %>% left_join(data_with_sp_cov)

}

#' Load grid from file path, crop extent and potentially mask
#'
#' @param file_path The full file path to a grid file
#' @param extent The extent of the grid file that should be loaded.
#' See \code{\link[raster]{crop}}.
#' @param shp If a (polygon) shp file provided, it will be used to
#' mask everything outside the specified polygons (can be slow,
#' depending on \code{shp} object. See \code{\link[raster]{mask}}.
#'
#' @importFrom raster raster crop mask
#' @importFrom sf as_Spatial
#'
#' @export
get_grid <- function(
  file_path,
  extent = extent(-180, 180, -90, 90),
  shp = NULL) {

  grid <- raster(file_path) %>% crop(extent)
  if (!is.null(shp)) {
    grid <- mask(grid, shp)
  }

  grid <- grid %>% as("SpatialGridDataFrame")
  grid@data[["year"]] <- str_extract_year(file_path)
  gc()

  grid

}


#' Import multiple grid data sets
#'
#' @param file_paths Full file paths to grid data.
#' @inheritParams get_grid
#'
#' @importFrom purrr map
#'
#' @export
get_grids <- function(
  file_paths,
  extent = extent(-180, 180, -90, 90),
  shp = NULL) {

  grids <- map(file_paths, ~ get_grid(.x, extent = extent, shp = shp))
  names(grids) <- get_file_name(file_paths)
  gc()

  grids

}
