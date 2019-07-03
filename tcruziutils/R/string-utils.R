#' Extract file name without without path and extension.
#'
#' @param file_path The full file path.
#'
#' @export
get_file_name <- function(file_path) {

  basename(tools::file_path_sans_ext(file_path))

}

#' Report years not available in grid files
#'
#' @param data Any data frame (whith year column \code{year_var}).
#' @param grid_files A vector of grid files. For each file the year will be
#' extracted using \code{str_extract_year}.
#' @param year_var The name of the variable stroring year information in
#' \code{data}.
#'
#' @importFrom dplyr pull
#' @importFrom mastergrids str_extract_year
#'
#' @export
missing_years_grid <- function(
  data,
  grid_files,
  year_var = "Start_year") {

  years_grid <- str_extract_year(grid_files) %>% sort()
  years_data <- data %>% pull(year_var) %>% unique() %>% na.omit() %>%
    as.integer() %>% sort()
  setdiff(years_data, years_grid)

}

#' List files by path and pattern
#'
#' A convenience wrapper around \code{list.files}
#'
#' @inheritParams base::list.files
#'
#' @export
list_files <- function(
  path,
  pattern,
  full.names = TRUE,
  recursive = TRUE) {

  list.files(path, pattern, full.names = full.names, recursive = recursive)

}
