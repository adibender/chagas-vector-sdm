#' Extract file name without without path and extension.
#'
#' @param file_path The full file path.
#' @importFrom tools file_path_sans_ext
#' @keywords internal
get_file_name <- function(file_path) {

  basename(file_path_sans_ext(file_path))

}

#' Extract year from character (vector)
#'
#' Given file path(s), extracts the year via regular expression and coverts
#' to integer.
#'
#' @param x A character (vector).
#' @importFrom stringr str_extract
#' @keywords internal
#' @export
str_extract_year <- function(x) {

  str_extract(basename(x), "20[0-9]+") %>% as.integer()

}


#' Define extent
#'
#' Depending on input, creates a \code{\link[raster]{extent}} object,
#' either from specified extent, the extent of a mask, or sets extent to the
#' maximum extent of a geographic map.
#'
#' @importFrom raster extent
#' @keywords internal
get_extent <- function(
  extent = NULL,
  mask   = NULL) {


  if (is.null(extent)) {
    if (!is.null(mask)) {
      extent <- extent(mask)
    }  else {
      extent <- extent(c(-180, 180, -90, 90))
    }
  }

  extent

}


#' Process names
#'
#' If provided, performs sanity checks, otherwise extracts names from
#' names of the raster layers provided in the \code{grid} argument.
#'
#' @importFrom purrr map_chr
#' @keywords internal
process_names <- function(names, files, data, grids) {

  assert_character(names, min.len = 1, null.ok = TRUE)

  if (is.null(names)) {
    names <- map_chr(grids, ~names(.x))
  } else {
    if (length(files) != length(names)) {
      stop("When specified, length of names should equal length of files.")
    }
    common_names <- intersect(names, colnames(data))
    if (!length(common_names) == 0) {
      warning(paste0("Column names ", paste0(common_names, collapse = ", "),
        "will be overwritten!"))
    }
  }

  names

}


#' Replace raster values by a column of a data frame
#'
#' Replaces values of \code{grid} with specified column at the coordinates of
#' raster that correspond to coordinates in \code{data} (column containing
#' coordinates are specified via \code{spat_form}). If \code{outside_na = TRUE}
#' (the default), \code{grid} values at coordinates outside of coordinates
#' in \code{data} are set to \code{NA}.
#'
#' @inheritParams get_grid_data
#' @param grid A Raster object or file to a GeoTIFF file that can be read in
#' via \code{raster}.
#' @param ... Further arguments passed to \code{read_grid} if \code{grid} is
#' a file instead of a raster object.
#' @importFrom checkmate checkFile
#' @importFrom raster cellFromXY values
#' @export
#' @keywords internal
df_to_grid <- function(
  data,
  grid,
  column,
  ...,
  spat_form = ~ longitude + latitude,
  reset_outside = TRUE,
  outside_value = NA) {

  if(testFile(grid, access = "r", extension = "tif")) {
    grid <- read_grid(grid, ...)
  }

  cells <- cellFromXY(grid, data[, all.vars(spat_form)])
  raster::values(grid)[cells] <- data[[column]]
  if(reset_outside) {
    raster::values(grid)[-cells] <- outside_value
  }

  grid

}
