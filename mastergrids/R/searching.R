#' List mastergrid paths
#'
#' A data frame containing names of mastergrid paths available in the
#' \code{mastergrid_paths} object and the respective values.
#'
#' @importFrom utils data
#' @export
list_paths <- function() {

  cbind.data.frame(
    name = names(mastergrids::mastergrid_paths),
    path = as.character(mastergrids::mastergrid_paths))
}

#' List all files in directory
#'
#'
#' @param ... Comma separated character strings that define parts of the file
# name. Only files within \code{path} will be kept, that contain all of the
#' parts defined in \code{...}. Specifications can also contain regular expressions.
#' @param file_extension Only files with this extension will be returned.
#' Usually we are only interested in GeoTiff files, thus \code{".tif$"} is
#' set by default.
#' @inheritParams base::list.files
#' @import checkmate
#' @importFrom purrr discard flatten map reduce
#' @examples
#' library(mastergrids)
#' # get mastergrids accessibility data (friction 2015)
#' search_files("friction", path = mastergrid_paths$accessibility_weiss)
#' search_files("Landcover", "Class02", "2005", "percentage",
#'   path = mastergrid_paths$landcover)
#' # search the entire mastergrid folder:
#' system.time({
#'   search_files("LST_Night", "2002", "Annual.max.5km")
#' })
#' @export
search_files <- function(
  ...,
  path           = "/media/z/mastergrids/",
  full.names     = TRUE,
  recursive      = TRUE,
  file_extension = ".tif$") {

  assert_directory(path, access = "r")
  assert_flag(full.names)
  assert_flag(recursive)
  assert_string(file_extension)

  dots <- list(...) %>% flatten()

  files <- list.files(
    path       = path,
    pattern    = file_extension,
    full.names = full.names,
    recursive  = recursive)
  ind_keep <- map(
    dots,
    ~ grep(.x, files, ignore.case = TRUE),
    .progress = .progress) %>%
    discard(~length(.x) == 0) %>%
    reduce(intersect)

  files[ind_keep]

}
