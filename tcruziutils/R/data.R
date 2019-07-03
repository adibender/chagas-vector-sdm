#' Spatial extent of Tcruzi analysis
#'
#' Can be used to create boundary box or used to clip spatial objects
"extent_tcruzi"

#' Boundary box for Tcruzi analysis
#'
#' Can be used to create boundary box or used to clip spatial objects
"bbox_tcruzi"

#' Spatial projection used for Tcruzi analysis
"proj_tcruzi"


#' Polygon plot of the endemic zone for the T.cruzi analysis
#'
#' Used as first layer in various thematic maps within the project
"p_endemic_zone"


#' Predifined color ramps
#'
#' Can be used to define color ramps for spatial plotting of environmental
#' covariates, e.g. as \code{palette} argument to \code{tmap::tm_*} objects.
#' @seealso \code{\link[grDevices]{colorRampPalette}}
#' @param n Number of colors returned.
#' @examples
#' WhRd(10)
"WhRd"

#' @describeIn WhRd Black Yellow color ramp for nighttime lights
"BlYl"
#' @describeIn WhRd Black White color ramp for accessibility
"BlWh"
#' @describeIn WhRd White Green color ramp for forests
"WhGr"
#' @describeIn WhRd White Brown color ramp
"WhBr"
#' @describeIn WhRd White Blue color ramp for rainfall
"WhBlue"


#' List of paths containing relevant environmental variables
#'
"tcruzi_paths"

#' List of regular expressions for specific GeoTIFF files
#'
"tcruzi_patterns"


#' A list of file paths to environmental variables
"tcruzi_queries"
