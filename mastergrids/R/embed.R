#' Extract coordinates from data and turn into SpatialPointsDataFrame
#'
#'
#' @importFrom raster projection
#' @importFrom sp coordinates 'coordinates<-'
#' @keywords internal
get_coords <- function(data, spat_form, proj) {

  spat_vars <- all.vars(spat_form)
  coords <- data[, c(spat_vars, ".temp_id")]
  coords <- na.omit(coords)
  coordinates(coords) <- spat_form
  projection(coords) <- proj

  coords

}


#' Extract grid values based on coordinates from data frame
#'
#' Given a data set that stores coordinates, extracts grid values from grid
#' data specified by the \code{files} argument. \code{get_grid_data_by_year}
#' does the same, but is based on multiple search queries and uses the
#' grid file closest to the respective year in \code{data}.
#'
#' @param data The data set that should be embeded by grid data. The
#' grid data is specified via \code{file = } argument. See \code{read_grid}
#' for details.
#' @inheritParams read_grids
#' @inheritParams raster::extract
#' @param ... Further arguments passed to \code{read_grid}. Note: usually
#' it is not necessary to provide an \code{extent} or a \code{mask}. In
#' fact it should not be done for \code{get_grid_data}, as it will slow down
#' the operation, as \code{\link[raster]{extract}} already only imports
#' grid data at coordinates specified in \code{data}.
#' @param spat_form A \code{formula} object. Specifies the longitude and
#' latitude column names in \code{data}. These will be used in combination with
#' \code{link[raster]{extract}} to embed \code{data} with grid data efficiently.
#' @param names A character string. If specified, will be used column name
#' that will store the grid information. Otherwise the default name will
#' be used.
#' @import checkmate
#' @importFrom raster extract projection 'projection<-'
#' @importFrom dplyr select right_join one_of everything
#' @importFrom purrr map map_dfc
#' @examples
#' library(mastergrids)
#' # create fake data
#' data <- sp::spsample(brazil, 6, "random") %>% as.data.frame()
#' # add grid data
#' file <- search_files("annual", "5km", "2002", "Class02", "perc",
#'   path=mastergrid_paths$landcover)
#' data %>% get_grid_data(file, spat_form = ~ x + y)
#' data %>% get_grid_data(file, spat_form =~ x + y, names = "Class00")
#'
#' # Also works correctly if some coordinates have missing values
#' data[1,] <- c(NA, -7.25)
#' data
#' data %>% get_grid_data(file, spat_form = ~x + y, names = "Class00")
#'
#' # Also works for multiple files
#' files <- search_files("annual", "5km", "2002", "Class02",
#'   path=mastergrid_paths$landcover)
#' data %>% get_grid_data(files, spat_form = ~x + y)
#' data %>% get_grid_data(files, spat_form = ~x + y, names = c("Adj", "Perc"))
#'
#' # Use get_grid_data_by_year to extract grid data based on year variable
#' queries <- list(
#'  search_files("acces.*weiss"),
#'  search_files("lst", "day", "annual", "5km", "mean.*mean"))
#' data$year <- 2000
#' data %>% get_grid_data_by_year(queries, spat_form = ~x + y,
#'  names = c("access", "lst_day"))
#' data$year <- 1997:2002
#' data %>% get_grid_data_by_year(queries, spat_form = ~x + y,
#'  names = c("access", "lst_day"))
#' # -> 1997 - 2000 the same lst_day as before (no values before 2000, so 2000
# is used)
#' # -> 2001 - 2003 uses the current lst_day values
#' # -> access always the same, as only one file available
#'
#' @export
get_grid_data <- function(
  data,
  files,
  ...,
  spat_form = ~ longitude + latitude,
  names     = NULL) {

  assert_data_frame(data, all.missing = TRUE, min.rows = 1, min.cols = 2)
  assert_file(files, access = "r")
  assert_character(names, null.ok = TRUE)
  assert_class(spat_form, "formula")
  spat_vars <- all.vars(spat_form)
  assert_character(spat_vars, len = 2)
  assert_subset(spat_vars, colnames(data))

  # get grid file
  grids <- map(files, function(file) read_grid(file, ...))
  names <- process_names(names, files, data, grids)
  # get coordinates from data
  cnames_data <- colnames(data)
  data[[".temp_id"]] <- seq_len(nrow(data))
  coords     <- get_coords(data, spat_form, projection(grids[[1]]))
  # extract grid values at coordinates and add to original data
  grid_vals <- map_dfc(grids, ~extract(.x, coords))
  colnames(grid_vals) <- names

  coords %>% as.data.frame() %>%
    cbind(grid_vals) %>%
    select(one_of(c(".temp_id", names))) %>%
    right_join(data, by = ".temp_id") %>%
    select(one_of(c(cnames_data, names)))

}

#' @rdname get_grid_data
#' @inherit get_grid_data
#' @param queries a list of search queries, each a results of \code{search_files}.
#' @param year_var The name of column in \code{data} that stores the year of
#' the observation.
#' @importFrom purrr map_dfr
#' @importFrom furrr future_map_dfr
#' @importFrom dplyr left_join select .data
#' @importFrom stats na.omit
#' @export
get_grid_data_by_year <- function(
  data,
  queries,
  year_var  = "year",
  ...,
  spat_form = ~ longitude + latitude,
  names     = NULL,
  .progress = interactive()) {

  assert_data_frame(data, all.missing = FALSE, min.rows = 1, min.cols = 3)
  assert_list(queries, types = "character")
  assert_character(names, null.ok = TRUE)
  assert_class(spat_form, "formula")
  spat_vars <- all.vars(spat_form)
  assert_character(spat_vars, len = 2)
  assert_subset(c(spat_vars, year_var), colnames(data))

  files_by_year <- group_files_by_year(queries, years = unique(data[[year_var]]))

  cnames <- colnames(data)
  data[[".temp_id2"]] <- seq_len(nrow(data))
  data_with_grid <- data[, c(spat_vars, ".temp_id2", year_var)] %>%
    na.omit() %>%
    split(f = .[[year_var]]) %>%
    future_map_dfr(
      ~get_grid_data(
        files = subset(files_by_year, year == .x[[year_var]][1])[["files"]][[1]],
        ...,
        spat_form = spat_form,
        names = names),
      .progress = .progress)
  data_with_grid %>% select(-one_of(c(spat_vars, year_var))) %>%
      right_join(data, by = ".temp_id2") %>%
      select(one_of(cnames), everything(), -.data$.temp_id2)

}


#' Find the file that matches the year most closely
#'
#'
#' @keywords internal
match_file_to_year <- function(files, year) {

  years_files   <- sort(str_extract_year(files))
  # catch case of synoptic files (without year specification)
  if (!length(years_files) == 0) {
    dist_year     <- abs(years_files - year)
  } else {
    dist_year <- 0
  }
  matching_file <- files[dist_year == min(dist_year)][1]

  matching_file

}

#' Group files by year
#'
#' See \code{match_file_by year} for details.
#' @importFrom tidyr nest unnest
#' @importFrom dplyr mutate select .data
#' @examples
#' queries <- list(
#'  search_files("acces.*weiss"),
#'  search_files("class00", "adj"))
#' mastergrids:::group_files_by_year(queries, years = 1997:2003)
#' @keywords internal
group_files_by_year <- function(
  queries,
  years,
  ...) {

  years <- as.numeric(years, any.missing = FALSE)
  assert_numeric(years)

  data.frame(year = as.numeric(years)) %>%
    mutate(files = map(.data$year,
      function(.year) {
        map_chr(queries, ~match_file_to_year(.x, .year)) %>%
          unique()
      })) %>%
    mutate(file_names = get_file_name(.data$files)) %>%
    select(.data$year, .data$file_names, .data$files)

}
