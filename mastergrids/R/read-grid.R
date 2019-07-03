#' Read (cropped/masked) grid data
#'
#' Given the file name of a GeoTiff file (including path), reads in the
#' grid data (as \code{RasterLayer}) object. If an extent is specified,
#' the file will be automatically cropped. If the spatial extent defined,
#' this can be more efficient (w.r.t. to time and memory) than reading in the
#' file and cropping it, #' as only the grid data within the defined extent will
#' be loaded into memory. Essentially, this is a wrapper around the respective
#' functionality in the \code{raster} package.
#'
#' @param file Use the \code{search_files} function to
#' easily obtain GeoTiff file names based on some prespecified criteria.
#' @param extent An object of class \code{\link[raster]{extent}}.
#' The spatial # extent of interest. Only the part of the grid
#' data within \code{extent} will be read in. If \code{NULL}, extent will be
#' determined from \code{mask}. If both NULL, the whole grid file will be read
#' in. Used as second argument in \code{\link[raster]{crop}}.
#' @param mask An object used to define a mask. All values of the
#' (cropped) grid outside the \code{mask} will be set to \code{NA}.
#' See \code{\link[raster]{mask}} for details.
#' @param ... Currently not used
#' @importFrom raster crop extent mask raster
#' @examples
#' library(mastergrids)
#' # search Landcover Class 00 files from 2002 with annual 5km resolution
#' # returns 2 files, adjacency and percentage
#' files <- search_files("landcover", "annual", "5km", "2002", "Class00")
#' # read in first file (spatial extent some subregion in south america):
#' system.time({
#'    grid1 <- read_grid(files[1]) %>% as("SpatialGridDataFrame")
#' })
#' sp::plot(grid1)
#' system.time({
#'    grid1 <- read_grid(files[1], extent = c(-61, -35, -35, -1)) %>%
#'      as("SpatialGridDataFrame")
#' })
#' sp::plot(grid1)
#' # Example with mask
#' # The grid file is first cropped w.r.t. to the extent of poly (extent(poly)),
#' # then all grid values outside the defined polygons are set to NA
#' grid2 <- read_grid(files[2], mask = brazil)
#' # Water outside the landmass was set to NA
#' sp::plot(grid2)
#' @export
read_grid <- function(
  file,
  extent = NULL,
  mask   = NULL,
  inverse = FALSE,
  ...) {

  assert_file(file, access = "r")

  extent <- get_extent(extent, mask)
  grid   <- raster(file)
  grid   <- grid %>% crop(extent)

  if (!is.null(mask)) {
    grid <- grid %>% mask(mask)
  }

  # add year attribute (can be used for merging with data)
  attr(grid@data, "year") <- str_extract_year(file)

  grid

}


#' Import multiple grid data sets
#'
#' @param files Full file paths to grid data.
#' @inheritParams read_grid
#' @inheritParams furrr::future_map
#'
#' @importFrom purrr map
#' @importFrom furrr future_map
#'
#' @export
read_grids <- function(
  files,
  extent = NULL,
  mask   = NULL,
  ...,
  .progress = interactive()) {

  assert_file(files, access = "r")

  grids <- future_map(
    files,
    ~ read_grid(.x, extent = extent, mask = mask, ...),
    .progress = .progress)

  grids

}
